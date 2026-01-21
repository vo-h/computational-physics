# Compiler settings
CC = clang
CFLAGS = -Wall -Wextra -std=c11 -g -I/opt/homebrew/opt/gsl/include
LDFLAGS = -L/opt/homebrew/opt/gsl/lib -lm -lgsl -lgslcblas -lcurl

# Directories
HF_SRC = hf/src
HF_INC = hf/include
HF_BUILD = hf/build
HF_BIN = hf/hf.exe
TEST_DIR = tests
TEST_BUILD = tests/build

# Source files for main HF program
HF_SOURCES = $(wildcard $(HF_SRC)/*.c)
HF_OBJECTS = $(patsubst $(HF_SRC)/%.c,$(HF_BUILD)/%.o,$(HF_SOURCES))

# Source files for tests (excluding cli.c)
HF_TEST_SOURCES = $(HF_SRC)/atom.c \
                  $(HF_SRC)/basis_stog.c \
                  $(HF_SRC)/molecule.c \
                  $(HF_SRC)/tensor.c
HF_TEST_OBJECTS = $(patsubst $(HF_SRC)/%.c,$(TEST_BUILD)/%.o,$(HF_TEST_SOURCES))

# Test executables
TESTS = test_atom test_tensor test_basis test_integrals

# Default target
all: hfc

# Create build directories
$(HF_BUILD):
	mkdir -p $(HF_BUILD)

$(TEST_BUILD):
	mkdir -p $(TEST_BUILD)

# Compile HF source files for main program
$(HF_BUILD)/%.o: $(HF_SRC)/%.c | $(HF_BUILD)
	$(CC) $(CFLAGS) -I$(HF_INC) -c $< -o $@

# Link main HF program
hfc: $(HF_OBJECTS)
	$(CC) $(HF_OBJECTS) -o $(HF_BIN) $(LDFLAGS)

# Compile HF source files for tests
$(TEST_BUILD)/%.o: $(HF_SRC)/%.c | $(TEST_BUILD)
	$(CC) $(CFLAGS) -I$(HF_INC) -c $< -o $@

# Compile test files
$(TEST_BUILD)/test_%.o: $(TEST_DIR)/test_%.c | $(TEST_BUILD)
	$(CC) $(CFLAGS) -I$(HF_INC) -c $< -o $@

# Build test executables
test_atom: $(TEST_BUILD)/test_atom.o $(HF_TEST_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $(TEST_BUILD)/$@ $(LDFLAGS)

test_tensor: $(TEST_BUILD)/test_tensor.o $(TEST_BUILD)/tensor.o
	$(CC) $(CFLAGS) $^ -o $(TEST_BUILD)/$@ $(LDFLAGS)

test_basis: $(TEST_BUILD)/test_basis.o $(TEST_BUILD)/basis_stog.o
	$(CC) $(CFLAGS) $^ -o $(TEST_BUILD)/$@ $(LDFLAGS)

test_integrals: $(TEST_BUILD)/test_integrals.o $(HF_TEST_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $(TEST_BUILD)/$@ $(LDFLAGS)

# Run all tests
test: $(TESTS)
	@echo ""
	@echo "=========================================="
	@echo "Running all tests..."
	@echo "=========================================="
	@cd tests && ./build/test_atom
	@echo ""
	@cd tests && ./build/test_tensor
	@echo ""
	@cd tests && ./build/test_basis
	@echo ""
	@cd tests && ./build/test_integrals
	@echo ""
	@echo "=========================================="
	@echo "All tests completed!"
	@echo "=========================================="

# Run individual tests
test-atom: test_atom
	cd tests && ./build/test_atom

test-tensor: test_tensor
	cd tests && ./build/test_tensor

test-basis: test_basis
	cd tests && ./build/test_basis

test-integrals: test_integrals
	cd tests && ./build/test_integrals

# Clean targets
clean-hfc:
	rm -rf $(HF_BUILD)/*.o $(HF_BIN)

clean-tests:
	rm -rf $(TEST_BUILD)

clean: clean-hfc clean-tests

.PHONY: all hfc test test-atom test-tensor test-basis test-integrals clean clean-hfc clean-tests