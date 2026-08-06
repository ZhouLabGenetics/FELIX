# designed by Kai, implemented by codex

CXX ?= g++
PKG_CONFIG ?= pkg-config

BIN_DIR ?= bin
STATIC_BIN_DIR ?= bin-static
BUILD_DIR ?= build
OBJ_DIR := $(BUILD_DIR)/felixla-objects
STATIC_OBJ_DIR := $(BUILD_DIR)/felixla-static-objects

FELIX_TARGET := $(BIN_DIR)/felixla
STATIC_FELIX_TARGET := $(STATIC_BIN_DIR)/felixla

FELIX_SRC := src/felixla.cpp
PACK_SRC := src/flare_subset_to_tractor_hybrid.cpp
VCF_SRC := src/tractor_hybrid_to_vcf.cpp
CMP_SRC := src/compare_vcfs.cpp
DOSAGE_SRC := src/tractor_dosage_vcf_to_hybrid.cpp
RFMIX_MSP_SRC := src/rfmix_msp_to_tractor_hybrid.cpp
EXTRACT_SRC := src/tractor_hybrid_extract_region.cpp
ADMIX_SRC := src/calc_tractor_admixture.cpp
FELIX_QUERY_SRC := src/felixla_query.cpp

TOOL_OBJS := \
	$(OBJ_DIR)/flare_subset_to_tractor_hybrid.o \
	$(OBJ_DIR)/tractor_hybrid_to_vcf.o \
	$(OBJ_DIR)/compare_vcfs.o \
	$(OBJ_DIR)/tractor_dosage_vcf_to_hybrid.o \
	$(OBJ_DIR)/rfmix_msp_to_tractor_hybrid.o \
	$(OBJ_DIR)/tractor_hybrid_extract_region.o \
	$(OBJ_DIR)/calc_tractor_admixture.o \
	$(OBJ_DIR)/felixla_query.o

STATIC_TOOL_OBJS := \
	$(STATIC_OBJ_DIR)/flare_subset_to_tractor_hybrid.o \
	$(STATIC_OBJ_DIR)/tractor_hybrid_to_vcf.o \
	$(STATIC_OBJ_DIR)/compare_vcfs.o \
	$(STATIC_OBJ_DIR)/tractor_dosage_vcf_to_hybrid.o \
	$(STATIC_OBJ_DIR)/rfmix_msp_to_tractor_hybrid.o \
	$(STATIC_OBJ_DIR)/tractor_hybrid_extract_region.o \
	$(STATIC_OBJ_DIR)/calc_tractor_admixture.o \
	$(STATIC_OBJ_DIR)/felixla_query.o

LEGACY_TOOL_NAMES := \
	flare_subset_to_tractor_hybrid \
	tractor_hybrid_to_vcf \
	compare_vcfs \
	tractor_dosage_vcf_to_hybrid \
	rfmix_msp_to_tractor_hybrid \
	tractor_hybrid_extract_region \
	calc_tractor_admixture \
	felixla_query \
	shapeit_chunks_to_tractor_args

PKG_HTSLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags htslib 2>/dev/null)
PKG_HTSLIB_LIBS := $(shell $(PKG_CONFIG) --libs htslib 2>/dev/null)
PKG_HTSLIB_STATIC_LIBS_RAW := $(shell $(PKG_CONFIG) --libs --static htslib 2>/dev/null)
PKG_HTSLIB_LIBDIR := $(shell $(PKG_CONFIG) --variable=libdir htslib 2>/dev/null)
PKG_HTSLIB_STATIC_ARCHIVE := $(PKG_HTSLIB_LIBDIR)/libhts.a
PKG_LIBDEFLATE_STATIC_LIBS := $(shell $(PKG_CONFIG) --libs --static libdeflate 2>/dev/null)
PKG_LIBDEFLATE_LIBDIR := $(shell $(PKG_CONFIG) --variable=libdir libdeflate 2>/dev/null)
PKG_LIBDEFLATE_STATIC_ARCHIVE := $(PKG_LIBDEFLATE_LIBDIR)/libdeflate.a
ifneq ($(wildcard $(PKG_LIBDEFLATE_STATIC_ARCHIVE)),)
PKG_LIBDEFLATE_STATIC_LIBS := $(PKG_LIBDEFLATE_STATIC_ARCHIVE)
endif
PKG_LIBLZMA_STATIC_LIBS := $(shell $(PKG_CONFIG) --libs --static liblzma 2>/dev/null)
PKG_LIBLZMA_LIBDIR := $(shell $(PKG_CONFIG) --variable=libdir liblzma 2>/dev/null)
PKG_LIBLZMA_STATIC_ARCHIVE := $(PKG_LIBLZMA_LIBDIR)/liblzma.a
ifneq ($(wildcard $(PKG_LIBLZMA_STATIC_ARCHIVE)),)
PKG_LIBLZMA_STATIC_LIBS := $(PKG_LIBLZMA_STATIC_ARCHIVE) $(filter-out -llzma -pthread -lpthread,$(PKG_LIBLZMA_STATIC_LIBS))
endif
PKG_LIBCURL_STATIC_LIBS := $(filter-out -lz,$(shell $(PKG_CONFIG) --libs --static libcurl 2>/dev/null))
HTSLIB_STATIC_DEPS := $(filter-out -lhts,$(PKG_HTSLIB_STATIC_LIBS_RAW))
ifneq ($(strip $(PKG_LIBDEFLATE_STATIC_LIBS)),)
HTSLIB_STATIC_DEPS := $(filter-out -ldeflate,$(HTSLIB_STATIC_DEPS)) $(PKG_LIBDEFLATE_STATIC_LIBS)
endif
ifneq ($(strip $(PKG_LIBLZMA_STATIC_LIBS)),)
HTSLIB_STATIC_DEPS := $(filter-out -llzma,$(HTSLIB_STATIC_DEPS)) $(PKG_LIBLZMA_STATIC_LIBS)
endif
ifneq ($(strip $(PKG_LIBCURL_STATIC_LIBS)),)
HTSLIB_STATIC_DEPS := $(filter-out -lcurl,$(HTSLIB_STATIC_DEPS)) $(PKG_LIBCURL_STATIC_LIBS)
endif
PKG_HTSLIB_STATIC_LIBS :=
ifneq ($(wildcard $(PKG_HTSLIB_STATIC_ARCHIVE)),)
PKG_HTSLIB_STATIC_LIBS := $(PKG_HTSLIB_STATIC_ARCHIVE) $(HTSLIB_STATIC_DEPS)
else
PKG_HTSLIB_STATIC_LIBS := $(PKG_HTSLIB_STATIC_LIBS_RAW)
endif
BREW_HTSLIB_PREFIX := $(shell brew --prefix htslib 2>/dev/null)
BREW_HTSLIB_CFLAGS :=
BREW_HTSLIB_LIBS :=
BREW_HTSLIB_STATIC_LIBS :=
ifneq ($(wildcard $(BREW_HTSLIB_PREFIX)/include/htslib/hts.h),)
BREW_HTSLIB_CFLAGS := -I$(BREW_HTSLIB_PREFIX)/include
BREW_HTSLIB_LIBS := -L$(BREW_HTSLIB_PREFIX)/lib -lhts
endif
ifneq ($(wildcard $(BREW_HTSLIB_PREFIX)/lib/libhts.a),)
BREW_HTSLIB_STATIC_LIBS := $(BREW_HTSLIB_PREFIX)/lib/libhts.a $(HTSLIB_STATIC_DEPS)
endif

HTSLIB_CFLAGS ?= $(if $(PKG_HTSLIB_CFLAGS),$(PKG_HTSLIB_CFLAGS),$(BREW_HTSLIB_CFLAGS))
HTSLIB_LIBS ?= $(if $(PKG_HTSLIB_LIBS),$(PKG_HTSLIB_LIBS),$(BREW_HTSLIB_LIBS))
ifeq ($(strip $(HTSLIB_LIBS)),)
HTSLIB_LIBS := -lhts
endif
HTSLIB_STATIC_LIBS ?= $(if $(PKG_HTSLIB_STATIC_LIBS),$(PKG_HTSLIB_STATIC_LIBS),$(BREW_HTSLIB_STATIC_LIBS))

CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic
STATIC_FULLY ?= 0
STATIC_LDFLAGS ?=
ifeq ($(STATIC_FULLY),1)
STATIC_LDFLAGS += -static
endif
CPPFLAGS += $(HTSLIB_CFLAGS)
LDLIBS += $(HTSLIB_LIBS)

.PHONY: all static clean clean-legacy-bin clean-legacy-static check-deps check-static-deps test test-static test-intense

all: clean-legacy-bin check-deps $(FELIX_TARGET)

static: clean-legacy-static check-static-deps $(STATIC_FELIX_TARGET)

clean-legacy-bin:
	rm -f $(addprefix $(BIN_DIR)/,$(LEGACY_TOOL_NAMES))
	@if [ -d "$(BIN_DIR)" ]; then \
		find "$(BIN_DIR)" -mindepth 1 -maxdepth 1 ! -name felixla -exec rm -R -f {} +; \
	fi

clean-legacy-static:
	rm -f $(addprefix $(STATIC_BIN_DIR)/,$(LEGACY_TOOL_NAMES))
	@if [ -d "$(STATIC_BIN_DIR)" ]; then \
		find "$(STATIC_BIN_DIR)" -mindepth 1 -maxdepth 1 ! -name felixla -exec rm -R -f {} +; \
	fi

check-deps:
	@printf '#include <htslib/hts.h>\n#include <htslib/vcf.h>\n' | \
		$(CXX) $(CPPFLAGS) -x c++ -E - >/dev/null 2>&1 || \
		( echo "ERROR: htslib headers not found. Install htslib or set HTSLIB_CFLAGS/HTSLIB_LIBS." >&2; exit 1 )

check-static-deps: check-deps
	@if [ -z "$(strip $(HTSLIB_STATIC_LIBS))" ]; then \
		echo "ERROR: htslib static library not found. Install libhts.a or set HTSLIB_STATIC_LIBS." >&2; \
		exit 1; \
	fi

$(FELIX_TARGET): $(FELIX_SRC) $(TOOL_OBJS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FELIX_SRC) $(TOOL_OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

$(STATIC_FELIX_TARGET): $(FELIX_SRC) $(STATIC_TOOL_OBJS)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FELIX_SRC) $(STATIC_TOOL_OBJS) -o $@ $(STATIC_LDFLAGS) $(HTSLIB_STATIC_LIBS)

$(OBJ_DIR)/flare_subset_to_tractor_hybrid.o: $(PACK_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_pack_main -c $< -o $@

$(OBJ_DIR)/tractor_hybrid_to_vcf.o: $(VCF_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_to_vcf_main -c $< -o $@

$(OBJ_DIR)/compare_vcfs.o: $(CMP_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_compare_main -c $< -o $@

$(OBJ_DIR)/tractor_dosage_vcf_to_hybrid.o: $(DOSAGE_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_dosage_main -c $< -o $@

$(OBJ_DIR)/rfmix_msp_to_tractor_hybrid.o: $(RFMIX_MSP_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_rfmix_main -c $< -o $@

$(OBJ_DIR)/tractor_hybrid_extract_region.o: $(EXTRACT_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_extract_main -c $< -o $@

$(OBJ_DIR)/calc_tractor_admixture.o: $(ADMIX_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_admixture_main -c $< -o $@

$(OBJ_DIR)/felixla_query.o: $(FELIX_QUERY_SRC)
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_query_main -c $< -o $@

$(STATIC_OBJ_DIR)/flare_subset_to_tractor_hybrid.o: $(PACK_SRC)
	mkdir -p $(STATIC_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_pack_main -c $< -o $@

$(STATIC_OBJ_DIR)/tractor_hybrid_to_vcf.o: $(VCF_SRC)
	mkdir -p $(STATIC_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_to_vcf_main -c $< -o $@

$(STATIC_OBJ_DIR)/compare_vcfs.o: $(CMP_SRC)
	mkdir -p $(STATIC_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_compare_main -c $< -o $@

$(STATIC_OBJ_DIR)/tractor_dosage_vcf_to_hybrid.o: $(DOSAGE_SRC)
	mkdir -p $(STATIC_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_dosage_main -c $< -o $@

$(STATIC_OBJ_DIR)/rfmix_msp_to_tractor_hybrid.o: $(RFMIX_MSP_SRC)
	mkdir -p $(STATIC_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_rfmix_main -c $< -o $@

$(STATIC_OBJ_DIR)/tractor_hybrid_extract_region.o: $(EXTRACT_SRC)
	mkdir -p $(STATIC_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_extract_main -c $< -o $@

$(STATIC_OBJ_DIR)/calc_tractor_admixture.o: $(ADMIX_SRC)
	mkdir -p $(STATIC_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_admixture_main -c $< -o $@

$(STATIC_OBJ_DIR)/felixla_query.o: $(FELIX_QUERY_SRC)
	mkdir -p $(STATIC_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Dmain=felixla_query_main -c $< -o $@

test: all
	BIN_DIR="$(abspath $(BIN_DIR))" bash tests/run_tiny.sh

test-intense: all
	python3 tests/run_keep_extract_intense.py --bin-dir "$(abspath $(BIN_DIR))"

test-static: static
	BIN_DIR="$(abspath $(STATIC_BIN_DIR))" bash tests/run_tiny.sh

clean:
	rm -rf $(BIN_DIR) $(STATIC_BIN_DIR) $(BUILD_DIR)
