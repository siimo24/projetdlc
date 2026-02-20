# ======================================================================
#  RSA Key Reconstruction — Makefile
# ======================================================================

CC       = gcc
CFLAGS   = -Wall -Wextra -O2 -std=c11 -I$(INCDIR)
LDFLAGS  = -lgmp -lssl -lcrypto -lm

SRCDIR   = src
INCDIR   = include
BUILDDIR = build
TARGET   = rsa_recon

# Source files (add new modules here)
SRCS     = $(SRCDIR)/main.c    \
           $(SRCDIR)/common.c  \
           $(SRCDIR)/keygen.c  \
           $(SRCDIR)/decay_sim.c \
           $(SRCDIR)/init_phase.c \
           $(SRCDIR)/heninger_shacham.c \
           $(SRCDIR)/henecka_may.c

OBJS     = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.o, $(SRCS))

# ---- Build rules ----

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "[BUILD] $(TARGET) built successfully"

$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# ---- Convenience targets ----

# Generate 1024-bit key with cold boot scenario
run: $(TARGET)
	./$(TARGET) --bits 1024 --preset cold --time 120 -v

# Quick test with 512-bit key
test: $(TARGET)
	./$(TARGET) --bits 512 --preset room --time 5 -v

# Simple decay mode (flat probability)
simple: $(TARGET)
	./$(TARGET) --bits 1024 --simple --known 0.30 --error 0.0

# Error model test (for HMM10)
error: $(TARGET)
	./$(TARGET) --bits 512 --simple --known 1.0 --error 0.15 -v

# Frozen scenario (liquid nitrogen)
frozen: $(TARGET)
	./$(TARGET) --bits 1024 --preset frozen --time 3600 -v

clean:
	rm -rf $(BUILDDIR) $(TARGET) keys/

.PHONY: all clean run test simple error frozen
