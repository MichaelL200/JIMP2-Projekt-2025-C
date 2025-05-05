# Kompilator i jego flagi
CC = cc
CFLAGS = -Wall -Wextra -Wuninitialized -Iinclude -I/usr/local/include/arpack
LDFLAGS = -L/usr/local/lib -larpack -lopenblas -lgfortran -lm

# Debugger i jego flagi
VALGRIND = time valgrind --suppressions=libgomp.supp --leak-check=full --show-leak-kinds=all --track-origins=yes

# Katalogi
SRC_DIR = src
OBJ_DIR = obj

# Nazwa pliku wykonywalnego
EXEC = graphdivider

# Pliki źródłowe i obiektowe
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC))

# Cel domyślny
all: $(EXEC)

# Tworzenie pliku wykonywalnego
$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $(EXEC) $(LDFLAGS)

# Kompilacja plików źródłowych
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

# Tworzenie katalogu obj/
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Czyszczenie plików pośrednich i końcowych
clean:
	rm -rf $(OBJ_DIR) $(EXEC)
clear: clean

# Czyszczenie i zbudowanie projektu
rebuild: clear all

# Cele do uruchamiania programu z Valgrind dla różnych plików wejściowych
g0: $(EXEC)
	$(VALGRIND) ./$(EXEC) input/g0.csrrg
0: $(EXEC)
	$(VALGRIND) ./$(EXEC) input/graf.csrrg
1: $(EXEC)
	$(VALGRIND) ./$(EXEC) input/graf1.csrrg
2: $(EXEC)
	$(VALGRIND) ./$(EXEC) input/graf2.csrrg
3: $(EXEC)
	$(VALGRIND) ./$(EXEC) input/graf3.csrrg
4: $(EXEC)
	$(VALGRIND) ./$(EXEC) input/graf4.csrrg
5: $(EXEC)
	$(VALGRIND) ./$(EXEC) input/graf5.csrrg
6: $(EXEC)
	$(VALGRIND) ./$(EXEC) input/graf6.csrrg

.PHONY: all clean clear rebuild g0 0 1 2 3 4 5 6
