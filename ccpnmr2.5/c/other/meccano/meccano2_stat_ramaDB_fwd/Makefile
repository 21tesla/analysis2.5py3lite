include ../../../environment.txt

DEBUG=no
# Paths=======================================================================#
DSRC = ../src
DINC = ../inc
DOBJ = ../obj
DEXE = ../bin
DPRO = .
#précision du compilateur=====================================================#
CC = gcc
#précision des options de compilation=========================================#
ifeq ($(DEBUG),yes)
	CFLAGS = $(GSL_INCLUDE_FLAGS) -Wall -W -Wwrite-strings -Wstrict-prototypes -Wuninitialized -Wunreachable-code -O3 -g -fPIC
else
	CFLAGS = $(GSL_INCLUDE_FLAGS) -Wall -W -O3 -fPIC
endif
#précision des options de création de liens=====================================#
LDFLAGS = $(GSL_LIB_FLAGS) $(GSL_LIB) -lm 
#nom de l'exécutable==========================================================#
EXEC = meccano2_stat_ramaDB_fwd
#liste de fichier source======================================================#
#Listage automatique de fichiers sources======================================#
SRC = $(wildcard $(DSRC)/*.c)
EXECSRC = $(DPRO)/meccano2_stat_ramaDB_main.c
#liste de objet===============================================================#
OBJ = $(SRC:.c=.o)
EXECOBJ = $(EXECSRC:.c=.o)

DEPEND= makedepend $(CFLAGS)
 
#Convention===================================================================#
#$@     : nom de la cicle
#$^     : liste des dépendances
#$<     : nom de la première dépendance
#$*     : nom du fichier sans suffixe ou extension
#@      : rend la commande muette ;'ex:@(CC)': 'make -s' en commande
 
all: meccano_main meccano_objects

meccano_main: $(DEXE)/$(EXEC)

$(DEXE)/$(EXEC): $(DPRO)/$(EXEC).o $(DPRO)/$(EXECOBJ) $(OBJ)
	@$(CC) -o $@ $^ $(LDFLAGS)
ifeq ($(DEBUG),yes)
	@echo "Génération en mode debug"
else
	@echo "Génération en mode release"
endif

meccano_objects: $(DPRO)/$(EXEC).o $(OBJ)

%.o : %.c
	@$(CC) -o $@ -c $< $(CFLAGS)

# what are the source dependencies
depend: $(SRCS)
	$(DEPEND) $(CFLAGS) $(SRCS)


.PHONY: clean mrproper

clean:
	@rm -f ../*/*.o ../*/*~

mrproper: clean
	@rm -f $(DEXE)/$(EXEC)


