#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "utilities.cpp"

int main(int argc, char* argv[]){
    /*Récupérer nom du problème*/
    try{
        std::string prob_name(argv[1]);

        // Déclaration problème
        OsiClpSolverInterface OSI_SOLVER;
        OsiClpSolverInterface SP_OSI_SOLVER;

        OSI_SOLVER.readMps(prob_name.c_str());
        double prob_UB = OSI_SOLVER.getInfinity();
        double prob_LB = -OSI_SOLVER.getInfinity();

        const unsigned int NBR_VARS = OSI_SOLVER.getNumCols();
        const unsigned int NBR_ROWS = OSI_SOLVER.getNumRows();

        // Set objective function sense (1 for min (default), -1 for max,)
        OSI_SOLVER.setObjSense(1);
        // Set log level (will also set underlying solver's log level)
        OSI_SOLVER.setLogLevel(3);

        // Récuparation fonction objective
        const double *OBJECTIVE = OSI_SOLVER.getObjCoefficients();
            double *c = new double[2];
            for(unsigned i=0; i<2; i++)
                c[i] = OBJECTIVE[i];

            double *f = new double[2];
            for(unsigned i=0; i<2; i++)
                f[i] = OBJECTIVE[i];

        // Récuparation matrice contraintes
        const CoinPackedMatrix *MATRIX = OSI_SOLVER.getMatrixByRow();

        // Récuparation second membre
        const double *ROW_UB = OSI_SOLVER.getRowUpper();
        const double *ROW_LB = OSI_SOLVER.getRowLower();



        //------------------------------
        std::cout <<std::endl << "Initial problem" << std::endl;
        displayProblem(&OSI_SOLVER);
        std::cout << std::endl;

        //-- Récupérer matrices A, B et b
        CoinPackedMatrix At;
        At.setDimensions(0, NBR_ROWS);
        for(unsigned i=0; i<2; i++){
            CoinPackedVector row;
            for(unsigned j=0; j<NBR_ROWS; j++)
                row.insert(j, MATRIX->getCoefficient(j,i));
            At.appendRow(row);
        }

        std::vector<std::vector<double> > B, b;
        B.resize(NBR_ROWS);
        b.resize(NBR_ROWS);
        for(unsigned i=0; i<NBR_ROWS; i++){
            B[i].resize(2);
            b[i].resize(1);
        }

        for(unsigned i=0; i<NBR_ROWS; i++){
            for(unsigned j=0; j<2; j++)
                B[i][j] = MATRIX->getCoefficient(i,j+2);
            b[i][0] = ROW_LB[i];
        }

        // Calcul fonction objective SP

        // test: y initial (5,5)
        std::vector<std::vector<double> > y_hat;
        y_hat.resize(2);
        for(unsigned i=0; i<y_hat.size(); i++){
            y_hat[i].resize(1);
        }

        for(unsigned i=0; i<y_hat.size(); i++){
            for(unsigned j=0; j<y_hat[0].size(); j++)
                y_hat[i][j] = 5;
        }

        std::vector<std::vector<double> > lambda;
        productMatrix(&B, &y_hat, &lambda);

        std::vector<std::vector<double> > SP_yhat;
        substractMatrix(&b, &lambda, &SP_yhat);

        //-- Générer sous problème 1
        generateSubProblem(&SP_OSI_SOLVER, &At, &SP_yhat, c, NBR_VARS-2, NBR_ROWS);
        displayProblem(&SP_OSI_SOLVER);
        SP_OSI_SOLVER.initialSolve();

        //-- Recuperer la solution ...
        const double * SOLUTION = SP_OSI_SOLVER.getColSolution();
        
        CoinAbsFltEq Equal;
         
        std::cout << "     Solution (non zero values only):" << std::endl << std::endl;
        for(unsigned int j=0; j<NBR_ROWS; ++j)
        {
            const double value = SOLUTION[j]; 
            
            if(!Equal(value,0))
            {
                std::ostringstream oss; oss << j; 
                std::string var_name("x" + oss.str());
                
                printf("       %8s   =  %3.2f \n", var_name.c_str(), value);
            }
        }



        std::cout << std::endl;

        if(c != NULL){delete[] c; c = NULL;}
        if(f != NULL){delete[] f; f = NULL;}
        if(ROW_UB != NULL){delete[] ROW_UB; ROW_UB = NULL;}
        if(ROW_LB != NULL){delete[] ROW_LB; ROW_LB = NULL;} 
        if(OBJECTIVE != NULL){delete[] OBJECTIVE; OBJECTIVE = NULL;} 
    }
    catch(std::exception & e){
        std::cerr << "\nException : " << e.what() << std::endl;
        std::cerr << "No file name\n" << std::endl;
        return -1;
    }
    return 0;
}