#ifndef _MYMACRO_H_
#define _MYMACRO_H_

#define MALLOC(_p,_t,_n)\
    do{ \
        (_p) = malloc ( (_n) * sizeof( (_t) )); \
        if( (_p) == NULL ) \
        { \
            fprintf(stderr,"Allocation impossible dans le fichier :%s ligne : %s",__FILE__,__LINE__);\
            exit(EXIT_FAILURE); \
        } \
    }while(0)

#define CALLOC(_p,_t,_n) \
    do{ \
        (_p) = calloc ( (_n) , sizeof( (_t) ) ); \
        if( (_p) == NULL ) \
        { \
            fprintf(stderr,"Allocation impossible dans le fichier :%s ligne : %s",__FILE__,__LINE__);\
            exit(EXIT_FAILURE); \
        } \
    }while(0)

#define REALLOC(_p,_t,_n)\
    do{ \
        (_t) * temp; \
        temp = realloc ( (_p) , (_n) * sizeof( (_t) ) ); \
        if( temp == NULL ) \
        { \
            fprintf(stderr,"Allocation impossible dans le fichier :%s ligne : %s",__FILE__,__LINE__);\
            free( (_p) ); \
            exit(EXIT_FAILURE); \
        } \
        else \
        { \
            (_p) = temp; \
        } \
    }while(0)

#define FREE(_p)\
    do{ \
        free( (_p) ); \
        (_p) = NULL; \
    }while(0)


#define SET(vector, i, value) gsl_vector_set(vector, i, value)
#define SET_ALL(vector, value) gsl_vector_set_all(vector, value)
#define GET(vector, i) gsl_vector_get(vector, i)
#define ISET(vector, i, value) gsl_vector_int_set(vector, i, value)
#define IGET(vector, i) gsl_vector_int_get(vector, i)
#define MSET(matrix, i, j, value) gsl_matrix_set(matrix, i, j, value)
#define MGET(matrix, i, j) gsl_matrix_get(matrix, i, j)
#define IMSET(matrix, i, j, value) gsl_matrix_int_set(matrix, i, j, value)
#define IMGET(matrix, i, j) gsl_matrix_int_get(matrix, i, j)

#define ALEA(a,b) gsl_ran_flat(r, a, b)    // Random number (uniform distribution) in [a b)
#define IALEA(a,b) (a + gsl_rng_uniform_int(r, b + 1 - a))  // Integer random number in [a b]

#define CHI2(calc,mes,sigma) gsl_pow_2(((calc)-(mes))/(sigma))
#define CHI2_FLAT(calc,mes,sigma) (gsl_pow_2(2.9846) * (1.0 - exp(-gsl_pow_2((((calc)-(mes))/(sigma))/2.9846))))

#endif
