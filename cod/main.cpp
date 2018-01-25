#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "utilities_matrix.cpp"
#include "utilities.cpp"

int main(int argc, char* argv[]){
    /*Récupérer nom du problème*/
    try{
        std::string prob_name(argv[1]);

        // Déclaration problème
        OsiClpSolverInterface OSI_SOLVER;
        OsiClpSolverInterface SP_OSI_SOLVER;
        OsiClpSolverInterface PPR_OSI_SOLVER;

        OSI_SOLVER.readMps(prob_name.c_str());
        prob_UB = OSI_SOLVER.getInfinity();
        prob_LB = -OSI_SOLVER.getInfinity();

        const unsigned int NBR_VARS = OSI_SOLVER.getNumCols();
        const unsigned int NBR_ROWS = OSI_SOLVER.getNumRows();

        NBR_VARS_Y = NBR_VARS/2;
        NBR_VARS_X = NBR_VARS - NBR_VARS_Y;

        // Set objective function sense (1 for min (default), -1 for max,)
        OSI_SOLVER.setObjSense(1);
        // Set log level (will also set underlying solver's log level)
        OSI_SOLVER.setLogLevel(3);

        // Récuparation fonction objective
        const double *OBJECTIVE = OSI_SOLVER.getObjCoefficients();
            double *c = new double[NBR_VARS_X];
            for(unsigned i=0; i<NBR_VARS_X; i++)
                c[i] = OBJECTIVE[i];

            double *f = new double[NBR_VARS_Y];
            for(unsigned i=0; i<NBR_VARS_Y; i++)
                f[i] = OBJECTIVE[i+NBR_VARS_X];

        // Récuparation matrice contraintes
        const CoinPackedMatrix *MATRIX = OSI_SOLVER.getMatrixByRow();

        // Récuparation second membre
        const double *ROW_UB = OSI_SOLVER.getRowUpper();
        const double *ROW_LB = OSI_SOLVER.getRowLower();



        //------------------------------
        std::cout <<std::endl << "\033[4mInitial problem\033[0m" << std::endl;
        displayProblem(&OSI_SOLVER);
        std::cout << std::endl;

        //-- Récupérer matrices A, B et b
        CoinPackedMatrix At;
        At.setDimensions(0, NBR_ROWS);
        for(unsigned i=0; i<NBR_VARS_X; i++){
            CoinPackedVector row;
            for(unsigned j=0; j<NBR_ROWS; j++)
                row.insert(j, MATRIX->getCoefficient(j,i));
            At.appendRow(row);
        }

        std::vector<std::vector<double> > B, b;
        B.resize(NBR_ROWS);
        b.resize(NBR_ROWS);
        for(unsigned i=0; i<NBR_ROWS; i++){
            B[i].resize(NBR_VARS_Y);
            b[i].resize(1);
        }

        for(unsigned i=0; i<NBR_ROWS; i++){
            for(unsigned j=0; j<NBR_VARS_Y; j++)
                B[i][j] = MATRIX->getCoefficient(i,j+NBR_VARS_X);
            b[i][0] = ROW_LB[i];
        }

        //------------------------------

        // Calcul fonction objective SP
        // test: yi initial 5
        std::vector<std::vector<double> > y_hat;
        y_hat.resize(NBR_VARS_Y);
        for(unsigned i=0; i<y_hat.size(); i++){
            y_hat[i].resize(1);
        }

        for(unsigned i=0; i<y_hat.size(); i++){
            for(unsigned j=0; j<y_hat[0].size(); j++)
                y_hat[i][j] = 5;
        }

        //-- Générer sous problème 1
        generateSubProblem(&SP_OSI_SOLVER, &At, &B, &b, &y_hat, c);
        //-- Générer problème principal 1
        generateMasterProblem(&PPR_OSI_SOLVER, &SP_OSI_SOLVER, f, &B, &b, &y_hat);
        updateYhat(&PPR_OSI_SOLVER, &y_hat);
        std::cout << prob_LB << " <= Optimal <= " << prob_UB << std::endl;

        unsigned iteration = 1;

        // std::cout << std::endl << "\033[4m Iteration n°" << iteration << " :\033[0m" << std::endl;
        // generateSubProblem(&SP_OSI_SOLVER, &At, &B, &b, &y_hat, c);
        // updateMasterProblem(&PPR_OSI_SOLVER, &SP_OSI_SOLVER, f, &B, &b, &y_hat);
        // updateYhat(&PPR_OSI_SOLVER, &y_hat);
        // std::cout << prob_LB << " <= Optimal <= " << prob_UB << std::endl;

        while((fabs(prob_UB - prob_LB) > 0.001) && iteration < 50){
            std::cout << std::endl << "\033[4m Iteration n°" << iteration << " :\033[0m" << std::endl;
            generateSubProblem(&SP_OSI_SOLVER, &At, &B, &b, &y_hat, c);
            updateMasterProblem(&PPR_OSI_SOLVER, &SP_OSI_SOLVER, f, &B, &b, &y_hat);
            updateYhat(&PPR_OSI_SOLVER, &y_hat);
            std::cout << prob_LB << " <= Optimal <= " << prob_UB << std::endl;

            iteration++;
            std::cout << std::endl;
        }

        // Nettoyage
        if(c != NULL){delete[] c; c = NULL;}
        if(f != NULL){delete[] f; f = NULL;}
        if(ROW_UB != NULL){delete[] ROW_UB; ROW_UB = NULL;}
        if(ROW_LB != NULL){delete[] ROW_LB; ROW_LB = NULL;} 
    }
    catch(std::exception & e){
        std::cerr << "\nException : " << e.what() << std::endl;
        std::cerr << "No file name\n" << std::endl;
        return -1;
    }
    return 0;
}