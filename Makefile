# Kompilator i jego flagi
CC = cc
CFLAGS = -Wall -Wextra -Wuninitialized -g -Iinclude
LDFLAGS = -lm

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

.PHONY: all clean clear
