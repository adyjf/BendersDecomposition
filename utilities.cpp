#include <iostream>

#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"


int productMatrix(const std::vector<std::vector<double> > *A, 
                const std::vector<std::vector<double> > *B, 
                std::vector<std::vector<double> > *C){
    unsigned m = A->size();
    unsigned n = (*A)[0].size();
    unsigned p = B->size();
    unsigned q = (*B)[0].size();

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

int substractMatrix(const std::vector<std::vector<double> > *A, 
                const std::vector<std::vector<double> > *B, 
                std::vector<std::vector<double> > *C){
    unsigned m = A->size();
    unsigned n = (*A)[0].size();
    unsigned p = B->size();
    unsigned q = (*B)[0].size();

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

void displayProblem(OsiClpSolverInterface *OSI_SOLVER){
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
}

void generateSubProblem(OsiClpSolverInterface *OSI_SOLVER,
                CoinPackedMatrix *A,
                const std::vector<std::vector<double> > *SP_yhat,
                double * c,
                unsigned nbr_vars, 
                unsigned nbr_rows){

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

    double * OBJECTIVE = new double[nbr_rows];
    for(unsigned i=0; i<nbr_rows; i++)
        OBJECTIVE[i] = (*SP_yhat)[i][0];

    OSI_SOLVER->loadProblem(*A, COL_LB, COL_UB, OBJECTIVE, ROW_LB, c);
    OSI_SOLVER->setObjSense(-1.);
    
    if(COL_UB != NULL){delete[] COL_UB; COL_UB = NULL;}
    if(COL_LB != NULL){delete[] COL_LB; COL_LB = NULL;}
    if(ROW_LB != NULL){delete[] ROW_LB; ROW_LB = NULL;} 
    if(OBJECTIVE != NULL){delete[] OBJECTIVE; OBJECTIVE = NULL;} 
}