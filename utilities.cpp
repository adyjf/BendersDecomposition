#include <iostream>

#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"

double prob_UB;
double prob_LB;
unsigned int NBR_VARS_X, NBR_VARS_Y;


int productMatrix(const std::vector<std::vector<double> > *A, 
                const std::vector<std::vector<double> > *B, 
                std::vector<std::vector<double> > *C){
    unsigned m = (unsigned)A->size();
    unsigned n = (unsigned)(*A)[0].size();
    unsigned p = (unsigned)B->size();
    unsigned q = (unsigned)(*B)[0].size();
    // std::cout << "m,n" << m << "," << n << std::endl;
    // std::cout << "p,q" << p << "," << q << std::endl;

    if(n == p){
        C->resize(m);
        for(unsigned i=0; i<m; i++)
            (*C)[i].resize(q);

        for(unsigned i=0; i<m; i++){
            for(unsigned j=0; j<q; j++){
                double tmp = 0;
                for(unsigned k=0; k<n; k++)
                        tmp += (*A)[i][k] * (*B)[k][j];
                (*C)[i][j] = tmp;
            }
        }

        return 1;
    }
    else
        return -1;
}

int productMatrixScalar(std::vector<std::vector<double> > *A, 
                double scalar){
    unsigned m = (unsigned)A->size();
    //unsigned n = (*A)[0].size();
    for(unsigned i=0; i<m; i++){
        for(unsigned j=0; j<1; j++){
            double tmp = (*A)[i][j] * scalar;
            (*A)[i][j] = tmp;
        }
    }

    return 1;
}

int substractMatrix(const std::vector<std::vector<double> > *A, 
                const std::vector<std::vector<double> > *B, 
                std::vector<std::vector<double> > *C){
    unsigned m = (unsigned)A->size();
    unsigned n = (unsigned)(*A)[0].size();
    unsigned p = (unsigned)B->size();
    unsigned q = (unsigned)(*B)[0].size();

    if(m == p && n == q){
        C->resize(m);
        for(unsigned i=0; i<m; i++)
            (*C)[i].resize(n);

        for(unsigned i=0; i<m; i++){
            for(unsigned j=0; j<q; j++){
                (*C)[i][j] = (*A)[i][j] - (*B)[i][j];
            }
        }

        return 1;
    }
    else
        return -1;
}

int addMatrix(const std::vector<std::vector<double> > *A, 
                const std::vector<std::vector<double> > *B, 
                std::vector<std::vector<double> > *C){
    unsigned m = (unsigned)A->size();
    unsigned n = (unsigned)(*A)[0].size();
    unsigned p = (unsigned)B->size();
    unsigned q = (unsigned)(*B)[0].size();

    if(m == p && n == q){
        C->resize(m);
        for(unsigned i=0; i<m; i++)
            (*C)[i].resize(n);

        for(unsigned i=0; i<m; i++){
            for(unsigned j=0; j<q; j++){
                (*C)[i][j] = (*A)[i][j] + (*B)[i][j];
            }
        }

        return 1;
    }
    else
        return -1;
}

void displayProblem(OsiClpSolverInterface *OSI_SOLVER){
    std::cout << std::endl;
    if(OSI_SOLVER->getObjSense() == 1)
        std::cout << "\tMinimize\t";
    else
        std::cout << "\tMaximize\t";

    const unsigned int nbr_vars = OSI_SOLVER->getNumCols();
    const unsigned int nbr_rows = OSI_SOLVER->getNumRows();
    // Récuparation fonction objective
    const double *objective = OSI_SOLVER->getObjCoefficients();

    // Récuparation matrice contraintes
    const CoinPackedMatrix *MATRIX = OSI_SOLVER->getMatrixByRow();

    // Récuparation second membre
    const double *row_ub = OSI_SOLVER->getRowUpper();
    const double *row_lb = OSI_SOLVER->getRowLower();

    for(unsigned int i=0; i<nbr_vars; i++){
        std::cout << objective[i] << "x" << i;
        if(i!=nbr_vars-1)
            std::cout << " + ";
    }
    std::cout << std::endl;

    std::cout << "\tsuch that :" << std::endl;
    for(unsigned int i=0; i<nbr_rows; i++){
        std::cout << "\t\t" << row_lb[i] << " <=\t";
        for(unsigned int j=0; j<nbr_vars; j++)
            std::cout << MATRIX->getCoefficient(i,j) << "\t";
        std::cout << "<= " <<row_ub[i] << std::endl;
    }
    std::cout << std::endl;
}

void modifySubProblem(OsiClpSolverInterface *SP_OSI_SOLVER){
    const unsigned int nbr_vars = SP_OSI_SOLVER->getNumCols();

    const double *SP_OBJECTIVE = SP_OSI_SOLVER->getObjCoefficients();
    CoinPackedVector row;
    for(unsigned j=0; j<nbr_vars; j++){
        row.insert(j, SP_OBJECTIVE[j]);
    }

    double * SP_ROW_UB = new double;
    double * SP_ROW_LB = new double;

    *SP_ROW_UB = 0;
    *SP_ROW_LB = - SP_OSI_SOLVER->getInfinity();

    SP_OSI_SOLVER->addRow(row, *SP_ROW_LB, *SP_ROW_UB);
}

void generateSubProblem(OsiClpSolverInterface *OSI_SOLVER,
                CoinPackedMatrix *A,
                std::vector<std::vector<double> > *B,
                std::vector<std::vector<double> > *b,
                std::vector<std::vector<double> > *y_hat,
                double * c){
    // Reinitialize solver
    OSI_SOLVER->reset();
    unsigned nbr_vars = A->getNumRows();
    unsigned nbr_rows = A->getNumCols();

    double * COL_UB = new double[nbr_rows];
    double * COL_LB = new double[nbr_rows];
    for(unsigned i=0; i<nbr_rows; ++i)
    {
        COL_LB[i] = 0;
        COL_UB[i] = OSI_SOLVER->getInfinity();
    }

    double * ROW_LB = new double[nbr_vars];
    for(unsigned i=0; i<nbr_vars; ++i)
        ROW_LB[i] = - OSI_SOLVER->getInfinity();

    std::vector<std::vector<double> > lambda;
    productMatrix(B, y_hat, &lambda);
    std::vector<std::vector<double> > SP_yhat;
    substractMatrix(b, &lambda, &SP_yhat);

    double * OBJECTIVE = new double[nbr_rows];
    for(unsigned i=0; i<nbr_rows; i++)
        OBJECTIVE[i] = SP_yhat[i][0];

    OSI_SOLVER->loadProblem(*A, COL_LB, COL_UB, OBJECTIVE, ROW_LB, c);
    OSI_SOLVER->setObjSense(-1.);
    
    if(COL_UB != NULL){delete[] COL_UB; COL_UB = NULL;}
    if(COL_LB != NULL){delete[] COL_LB; COL_LB = NULL;}
    if(ROW_LB != NULL){delete[] ROW_LB; ROW_LB = NULL;} 
    if(OBJECTIVE != NULL){delete[] OBJECTIVE; OBJECTIVE = NULL;}
}

void generateMasterProblem(OsiClpSolverInterface *PPR_OSI_SOLVER,
                OsiClpSolverInterface *SP_OSI_SOLVER,
                double * f,
                std::vector<std::vector<double> > *B,
                std::vector<std::vector<double> > *b,
                std::vector<std::vector<double> > *y_hat){

    const unsigned int nbr_vars = SP_OSI_SOLVER->getNumCols();
    const unsigned int nbr_rows = SP_OSI_SOLVER->getNumRows();
    // std::cout << nbr_vars << "\t" << nbr_cols;
    displayProblem(SP_OSI_SOLVER);
    SP_OSI_SOLVER->initialSolve();

    //-- Recuperer la solution du sous probleme
    const double * sp_solution = SP_OSI_SOLVER->getColSolution();

    std::vector<std::vector<double> > sp_solution_;
    sp_solution_.resize(1);
    for(unsigned i=0; i<1; i++)
        sp_solution_[i].resize(nbr_vars);

    for(unsigned i=0; i<sp_solution_.size(); i++){
        for(unsigned j=0; j<sp_solution_[i].size(); j++)
            sp_solution_[i][j] = sp_solution[j];
    }

    CoinAbsFltEq Equal;
     
    std::cout << "     Solution (non zero values only):" << std::endl << std::endl;
    for(unsigned int j=0; j<nbr_vars; ++j){
        const double value = sp_solution[j]; 
        
        if(!Equal(value,0)){
            std::ostringstream oss; oss << j; 
            std::string var_name("x" + oss.str());
            
            printf("       %8s   =  %3.2f \n", var_name.c_str(), value);
        }
    }

    // Déclaration problème maitre
    // Fonction objective
    double * PPR_OBJECTIVE = new double[nbr_rows+1];
    for(unsigned j=0; j<nbr_rows; j++){
        PPR_OBJECTIVE[j] = 0;
    }
    PPR_OBJECTIVE[nbr_rows] = 1;

    CoinPackedMatrix PPR_MATRIX;
        PPR_MATRIX.setDimensions(0,nbr_rows+1);
    double * PPR_ROW_UB = new double;
    double * PPR_ROW_LB = new double;

    if(SP_OSI_SOLVER->isProvenOptimal()){
        // Maj prob_UB
        std::vector<std::vector<double> > f_;
        f_.resize(1);
        f_[0].resize(nbr_rows);
        for(unsigned j=0; j<nbr_rows; j++)
            f_[0][j] = f[j];
        std::vector<std::vector<double> > f_yhat;
        productMatrix(&f_, y_hat, &f_yhat);
        prob_UB = SP_OSI_SOLVER->getObjValue() + f_yhat[0][0];

        // Matrice des contraintes
        std::vector<std::vector<double> > lambda;
        productMatrix(&sp_solution_, B, &lambda);
        // productMatrixScalar(&lambda, -1);

        CoinPackedVector row;
        for(unsigned j=0; j<nbr_rows; j++){
            const double value = f[j] - lambda[0][j];
            row.insert(j, value);
        }
        row.insert(nbr_rows, -1);
        PPR_MATRIX.appendRow(row);

        // Second Membre
        std::vector<std::vector<double> > ppr_row_ub_;    
        productMatrix(&sp_solution_, b, &ppr_row_ub_);

        *PPR_ROW_UB = -ppr_row_ub_[0][0];

        *PPR_ROW_LB = - PPR_OSI_SOLVER->getInfinity();
    }
    else{
        std::cout << "Unbounded TODO" << std::endl;
        modifySubProblem(SP_OSI_SOLVER);
        SP_OSI_SOLVER->initialSolve();

        //-- Recuperer la solution du sous probleme modifié
        const double * sp_solution = SP_OSI_SOLVER->getColSolution();

        std::vector<std::vector<double> > sp_solution_;
        sp_solution_.resize(1);
        for(unsigned i=0; i<1; i++)
            sp_solution_[i].resize(nbr_vars);

        for(unsigned i=0; i<sp_solution_.size(); i++){
            for(unsigned j=0; j<sp_solution_[i].size(); j++)
                sp_solution_[i][j] = sp_solution[j];
        }

        // Matrice des contraintes
        std::vector<std::vector<double> > lambda;
        productMatrix(&sp_solution_, B, &lambda);
        // productMatrixScalar(&lambda, -1);

        CoinPackedVector row;
        for(unsigned j=0; j<nbr_rows; j++){
            const double value = lambda[0][j];
            row.insert(j, value);
        }
        row.insert(nbr_rows, -1);
        PPR_MATRIX.appendRow(row);

        // Second Membre
        std::vector<std::vector<double> > ppr_row_ub_;    
        productMatrix(&sp_solution_, b, &ppr_row_ub_);

        *PPR_ROW_UB = -ppr_row_ub_[0][0];

        *PPR_ROW_LB = - PPR_OSI_SOLVER->getInfinity();
    }

    // Borne des variables
    double * PPR_COL_UB = new double[nbr_rows+1];
    double * PPR_COL_LB = new double[nbr_rows+1];

    for(unsigned int i=0; i<nbr_rows; ++i){
        PPR_COL_LB[i] = 0;
        PPR_COL_UB[i] = PPR_OSI_SOLVER->getInfinity();
    }
    PPR_COL_LB[nbr_rows] = -PPR_OSI_SOLVER->getInfinity();
    PPR_COL_UB[nbr_rows] = PPR_OSI_SOLVER->getInfinity();

    PPR_OSI_SOLVER->loadProblem(PPR_MATRIX, PPR_COL_LB, PPR_COL_UB, PPR_OBJECTIVE, PPR_ROW_LB, PPR_ROW_UB);
    PPR_OSI_SOLVER->setObjSense(1);

    // Nettoyage
    if(PPR_COL_UB != NULL){delete[] PPR_COL_UB; PPR_COL_UB = NULL;}
    if(PPR_COL_LB != NULL){delete[] PPR_COL_LB; PPR_COL_LB = NULL;}
    if(PPR_ROW_UB != NULL){delete PPR_ROW_UB; PPR_ROW_UB = NULL;}
    if(PPR_ROW_LB != NULL){delete PPR_ROW_LB; PPR_ROW_LB = NULL;} 
    if(PPR_OBJECTIVE != NULL){delete[] PPR_OBJECTIVE; PPR_OBJECTIVE = NULL;}
}

void updateYhat(OsiClpSolverInterface *PPR_OSI_SOLVER,
                std::vector<std::vector<double> > *y_hat){
    displayProblem(PPR_OSI_SOLVER);
    PPR_OSI_SOLVER->initialSolve();

    //-- Recuperer la solution du probleme principal
    const double * ppr_solution = PPR_OSI_SOLVER->getColSolution();
    CoinAbsFltEq Equal;

    std::cout << "     Solution (non zero values only):" << std::endl << std::endl;
    for(unsigned j=0; j<(unsigned)PPR_OSI_SOLVER->getNumCols(); ++j){
        const double value = ppr_solution[j]; 

        if(!Equal(value,0)){
        std::ostringstream oss; oss << j; 
        std::string var_name("x" + oss.str());

        printf("       %8s   =  %3.2f \n", var_name.c_str(), value);
        }
    }
    // Mise à jour LB du probleme
    prob_LB = PPR_OSI_SOLVER->getObjValue();

    for(unsigned i=0; i<y_hat->size(); i++){
        for(unsigned j=0; j<(*y_hat)[0].size(); j++)
            (*y_hat)[i][j] = ppr_solution[i];
    }
}

void updateMasterProblem(OsiClpSolverInterface *PPR_OSI_SOLVER,
                OsiClpSolverInterface *SP_OSI_SOLVER,
                double * f,
                std::vector<std::vector<double> > *B,
                std::vector<std::vector<double> > *b,
                std::vector<std::vector<double> > *y_hat){

    const unsigned int nbr_vars = SP_OSI_SOLVER->getNumCols();
    const unsigned int nbr_rows = SP_OSI_SOLVER->getNumRows();
    // std::cout << nbr_vars << "\t" << nbr_cols;
    displayProblem(SP_OSI_SOLVER);
    SP_OSI_SOLVER->initialSolve();

    //-- Recuperer la solution du sous probleme
    const double * sp_solution = SP_OSI_SOLVER->getColSolution();

    std::vector<std::vector<double> > sp_solution_;
    sp_solution_.resize(1);
    for(unsigned i=0; i<1; i++)
        sp_solution_[i].resize(nbr_vars);

    for(unsigned i=0; i<sp_solution_.size(); i++){
        for(unsigned j=0; j<sp_solution_[i].size(); j++)
            sp_solution_[i][j] = sp_solution[j];
    }

    CoinAbsFltEq Equal;
     
    std::cout << "     Solution (non zero values only):" << std::endl << std::endl;
    for(unsigned int j=0; j<nbr_vars; ++j){
        const double value = sp_solution[j]; 
        
        if(!Equal(value,0)){
            std::ostringstream oss; oss << j; 
            std::string var_name("x" + oss.str());
            
            printf("       %8s   =  %3.2f \n", var_name.c_str(), value);
        }
    }

    // Déclaration problème maitre
    double * PPR_ROW_UB = new double;
    double * PPR_ROW_LB = new double;

    if(SP_OSI_SOLVER->isProvenOptimal()){
        // Maj prob_UB
        std::vector<std::vector<double> > f_;
        f_.resize(1);
        f_[0].resize(nbr_rows);
        for(unsigned j=0; j<nbr_rows; j++)
            f_[0][j] = f[j];
        std::vector<std::vector<double> > f_yhat;
        productMatrix(&f_, y_hat, &f_yhat);
        prob_UB = SP_OSI_SOLVER->getObjValue() + f_yhat[0][0];

        // Matrice des contraintes
        CoinPackedVector row;
        std::vector<std::vector<double> > lambda;
        productMatrix(&sp_solution_, B, &lambda);
        for(unsigned j=0; j<nbr_rows; j++){
            const double value = f[j] - lambda[0][j];
            row.insert(j, value);
        }
        row.insert(nbr_rows, -1);

        // Second Membre
        std::vector<std::vector<double> > ppr_row_ub_;    
        productMatrix(&sp_solution_, b, &ppr_row_ub_);

        *PPR_ROW_UB = -ppr_row_ub_[0][0];

        *PPR_ROW_LB = - PPR_OSI_SOLVER->getInfinity();

        PPR_OSI_SOLVER->addRow(row, *PPR_ROW_LB, *PPR_ROW_UB);
    }
    else{
        std::cout << "Unbounded TODO" << std::endl;
        modifySubProblem(SP_OSI_SOLVER);
        SP_OSI_SOLVER->initialSolve();

        //-- Recuperer la solution du sous probleme modifié
        const double * sp_solution = SP_OSI_SOLVER->getColSolution();

        std::vector<std::vector<double> > sp_solution_;
        sp_solution_.resize(1);
        for(unsigned i=0; i<1; i++)
            sp_solution_[i].resize(nbr_vars);

        for(unsigned i=0; i<sp_solution_.size(); i++){
            for(unsigned j=0; j<sp_solution_[i].size(); j++)
                sp_solution_[i][j] = sp_solution[j];
        }

        // Matrice des contraintes
        std::vector<std::vector<double> > lambda;
        productMatrix(&sp_solution_, B, &lambda);
        // productMatrixScalar(&lambda, -1);

        CoinPackedVector row;
        for(unsigned j=0; j<nbr_rows; j++){
            const double value = lambda[0][j];
            row.insert(j, value);
        }
        row.insert(nbr_rows, -1);

        // Second Membre
        std::vector<std::vector<double> > ppr_row_ub_;    
        productMatrix(&sp_solution_, b, &ppr_row_ub_);

        *PPR_ROW_UB = -ppr_row_ub_[0][0];

        *PPR_ROW_LB = - PPR_OSI_SOLVER->getInfinity();

        PPR_OSI_SOLVER->addRow(row, *PPR_ROW_LB, *PPR_ROW_UB);
    }
    
    // Nettoyage
    if(PPR_ROW_UB != NULL){delete PPR_ROW_UB; PPR_ROW_UB = NULL;}
    if(PPR_ROW_LB != NULL){delete PPR_ROW_LB; PPR_ROW_LB = NULL;}
}