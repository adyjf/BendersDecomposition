// productMatrix : fait le produit de la matrice A avec la matrice B et met le résultat dans C
int productMatrix(const std::vector<std::vector<double> > *A, 
                const std::vector<std::vector<double> > *B, 
                std::vector<std::vector<double> > *C){
    unsigned m = (unsigned)A->size();
    unsigned n = (unsigned)(*A)[0].size();
    unsigned p = (unsigned)B->size();
    unsigned q = (unsigned)(*B)[0].size();

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

// substractMatrix : soustraction de matrices
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

// addMatrix : addition de matrices
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