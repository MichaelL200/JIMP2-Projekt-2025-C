# Kompilator i jego flagi
CC = cc
CFLAGS = -Wall -Wextra -Wuninitialized -g -Iinclude
LDFLAGS = -lm

# Katalogi
SRC_DIR = src
OBJ_DIR = obj

# Nazwa pliku wykonywalnego
EXEC = graphdivider
LINK = gd

# Pliki źródłowe i obiektowe
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC))

# Cel domyślny
all: $(EXEC) $(LINK)

# Tworzenie pliku wykonywalnego
$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $(EXEC) $(LDFLAGS)

# Tworzenie dowiązania symbolicznego
$(LINK): $(EXEC)
	ln -sf $(EXEC) $(LINK)

# Kompilacja plików źródłowych
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

# Tworzenie katalogu obj/
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Czyszczenie plików pośrednich i końcowych
clean:
	rm -rf $(OBJ_DIR) $(EXEC) $(LINK)
clear: clean

.PHONY: all clean clear
