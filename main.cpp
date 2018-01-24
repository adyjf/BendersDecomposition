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
        OsiClpSolverInterface PPR_OSI_SOLVER;

        OSI_SOLVER.readMps(prob_name.c_str());
        prob_UB = OSI_SOLVER.getInfinity();
        prob_LB = -OSI_SOLVER.getInfinity();

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
                f[i] = OBJECTIVE[i+2];

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
        //------------------------------

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

        //-- Générer sous problème 1
        generateSubProblem(&SP_OSI_SOLVER, &At, &B, &b, &y_hat, c);
        generateMasterProblem(&PPR_OSI_SOLVER, &SP_OSI_SOLVER, f, &B, &b);
        updateYhat(&PPR_OSI_SOLVER, &y_hat);

        //-- Générer sous problème 2
        generateSubProblem(&SP_OSI_SOLVER, &At, &B, &b, &y_hat, c);
        updateMasterProblem(&PPR_OSI_SOLVER, &SP_OSI_SOLVER, f, &B, &b);
        updateYhat(&PPR_OSI_SOLVER, &y_hat);

        //-- Générer sous problème 3
        generateSubProblem(&SP_OSI_SOLVER, &At, &B, &b, &y_hat, c);
        updateMasterProblem(&PPR_OSI_SOLVER, &SP_OSI_SOLVER, f, &B, &b);
        updateYhat(&PPR_OSI_SOLVER, &y_hat);

        //-- Générer sous problème 4
        generateSubProblem(&SP_OSI_SOLVER, &At, &B, &b, &y_hat, c);
        updateMasterProblem(&PPR_OSI_SOLVER, &SP_OSI_SOLVER, f, &B, &b);
        updateYhat(&PPR_OSI_SOLVER, &y_hat);

        //-- Générer sous problème 5
        generateSubProblem(&SP_OSI_SOLVER, &At, &B, &b, &y_hat, c);
        updateMasterProblem(&PPR_OSI_SOLVER, &SP_OSI_SOLVER, f, &B, &b);
        updateYhat(&PPR_OSI_SOLVER, &y_hat);

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