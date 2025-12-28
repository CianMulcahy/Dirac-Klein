# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Iinclude

# Source and object files
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)

# Executable name
TARGET = Dirac

# Default rule: build executable
all: $(TARGET)

# Link objects to create executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET)

# Compile each .cpp into .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f src/*.o $(TARGET)