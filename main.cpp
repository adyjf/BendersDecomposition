#include <iostream>
#include <string>
#include <vector>

#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"

int main(int argc, char* argv[]){
    /*Récupérer nom du problème*/
    try{
        std::string prob_name(argv[1]);

        OsiClpSolverInterface OSI_SOLVER;
        CbcModel CBCMODEL(OSI_SOLVER);

        CBCMODEL.solver()->readMps(prob_name.c_str());
        const unsigned int NBR_VARS = CBCMODEL.solver()->getNumCols();
        const unsigned int NBR_ROWS = CBCMODEL.solver()->getNumRows();

        // à update: 
        // Set objective function sense (1 for min (default), -1 for max,)
        // CBCMODEL.solver()->setObjSense(1);
        // Set log level (will also set underlying solver's log level)
        // CBCMODEL.setLogLevel(3);
        std::cout << "ok" << std::endl;
        const double * z = CBCMODEL.solver()->getRowUpper();

        std::cout << "ok1" << std::endl;
        for(unsigned int i=0; i<NBR_ROWS; i++){
            std::cout << z[i] << std::endl;
        }
        std::cout << "ok2" << std::endl;

        const CoinPackedMatrix *MATRIX = CBCMODEL.solver()->getMatrixByRow();
        std::vector<std::vector<double> > x;
        x.resize(NBR_ROWS);
        for(unsigned int i=0; i<NBR_ROWS; i++){
            x[i].resize(NBR_VARS);
        }

        for(unsigned int i=0; i<NBR_ROWS; i++){
            for(unsigned int j=0; j<NBR_VARS; j++){
                x[i][j] = MATRIX->getCoefficient(i,j);
                std::cout << x[i][j] << "\t";
            }
            std::cout << std::endl;
        }
    }
    catch(std::exception & e){
        std::cerr << "\nException : " << e.what() << std::endl;
        std::cerr << "No file name\n" << std::endl;
        return -1;
    }
    return 0;
}