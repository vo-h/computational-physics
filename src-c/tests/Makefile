# Compiler settings
CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -g -I../hf/include -I/opt/homebrew/opt/gsl/include
LDFLAGS = -L/opt/homebrew/opt/gsl/lib -lm -lgsl -lgslcblas -lcurl

# Directories
HF_SRC_DIR = ../hf/src
TEST_DIR = .
BUILD_DIR = build

# Source files from hf directory
HF_SOURCES = $(HF_SRC_DIR)/atom.c \
             $(HF_SRC_DIR)/basis_stog.c \
             $(HF_SRC_DIR)/molecule.c \
             $(HF_SRC_DIR)/tensor.c

HF_OBJECTS = $(patsubst $(HF_SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(HF_SOURCES))

# Test executables
TESTS = test_atom test_tensor test_basis test_integrals

# Default target
all: $(BUILD_DIR) $(TESTS)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile HF source files
$(BUILD_DIR)/%.o: $(HF_SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Build test executables
test_atom: $(BUILD_DIR)/test_atom.o $(HF_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $(BUILD_DIR)/$@ $(LDFLAGS)

test_tensor: $(BUILD_DIR)/test_tensor.o $(BUILD_DIR)/tensor.o
	$(CC) $(CFLAGS) $^ -o $(BUILD_DIR)/$@ $(LDFLAGS)

test_basis: $(BUILD_DIR)/test_basis.o $(BUILD_DIR)/basis_stog.o
	$(CC) $(CFLAGS) $^ -o $(BUILD_DIR)/$@ $(LDFLAGS)

test_integrals: $(BUILD_DIR)/test_integrals.o $(HF_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $(BUILD_DIR)/$@ $(LDFLAGS)

# Compile test files
$(BUILD_DIR)/test_%.o: $(TEST_DIR)/test_%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Run all tests
test: all
	@echo ""
	@echo "=========================================="
	@echo "Running all tests..."
	@echo "=========================================="
	@./$(BUILD_DIR)/test_atom
	@echo ""
	@./$(BUILD_DIR)/test_tensor
	@echo ""
	@./$(BUILD_DIR)/test_basis
	@echo ""
	@./$(BUILD_DIR)/test_integrals
	@echo ""
	@echo "=========================================="
	@echo "All tests completed!"
	@echo "=========================================="

# Run individual tests
test-atom: test_atom
	./$(BUILD_DIR)/test_atom

test-tensor: test_tensor
	./$(BUILD_DIR)/test_tensor

test-basis: test_basis
	./$(BUILD_DIR)/test_basis

test-integrals: test_integrals
	./$(BUILD_DIR)/test_integrals

# Clean up
clean:
	rm -rf $(BUILD_DIR)

# Phony targets
.PHONY: all test test-atom test-tensor test-basis test-integrals clean
