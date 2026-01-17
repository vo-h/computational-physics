CC=clang
# CFLAGS=-Wall -Wextra
CFLAGS=
HF_SRC=src-c/hf/src
HF_INC=src-c/hf/include
HF_BUILD=src-c/hf/build
HF_BIN=src-c/hf/hf.exe

# Source files
HF_SOURCES=$(wildcard $(HF_SRC)/*.c)
HF_OBJECTS=$(patsubst $(HF_SRC)/%.c,$(HF_BUILD)/%.o,$(HF_SOURCES))

# Create build directory if it doesn't exist
$(HF_BUILD):
	mkdir -p $(HF_BUILD)

# Compile source files to object files in build directory
$(HF_BUILD)/%.o: $(HF_SRC)/%.c | $(HF_BUILD)
	$(CC) $(CFLAGS) -I$(HF_INC) -I/opt/homebrew/Cellar/gsl/2.8/include -c $< -o $@

# Link all object files to create final program
hfc: $(HF_OBJECTS)
	$(CC) $(HF_OBJECTS) -o $(HF_BIN) -lm -lcurl -I/opt/homebrew/Cellar/gsl/2.8/include -L/opt/homebrew/Cellar/gsl/2.8/lib -lgsl -lgslcblas

.PHONY: hfc clean-hfc

clean-hfc:
	rm -rf $(HF_BUILD)/*.o $(HF_BIN)