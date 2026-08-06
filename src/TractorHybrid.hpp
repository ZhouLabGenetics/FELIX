// TractorHybrid.hpp
// designed by Kai, implemented by codex
//
// FELIX step2 reader for the ancestry-aware packed hybrid format.

#ifndef TRACTOR_HYBRID_HPP
#define TRACTOR_HYBRID_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace TRACTOR_HYBRID {

struct RareCarrierPacked {
    uint32_t pos_index;
    uint32_t anc_hap;
};

class TractorHybridClass {
public:
    struct Meta;
    struct VariantRecord;
    class MksStream;

    TractorHybridClass(
        std::string prefix,
        std::vector<std::string> sampleInModel
    );

    ~TractorHybridClass();

    void getOneMarker(
        uint64_t& t_gIndex_prev,
        uint64_t& t_gIndex,
        std::string& t_ref,
        std::string& t_alt,
        std::string& t_marker,
        uint32_t& t_pd,
        std::string& t_chr,
        double& t_altFreq,
        double& t_altCounts,
        double& t_missingRate,
        double& t_imputeInfo,
        bool& t_isOutputIndexForMissing,
        std::vector<uint>& t_indexForMissing,
        bool& t_isOnlyOutputNonZero,
        std::vector<uint>& t_indexForNonZero,
        bool& t_isBoolRead,
        arma::vec& t_GVec,
        bool t_isImputation
    );

    void getOneMarkerAncestry(
        uint64_t& t_gIndex_prev,
        uint64_t& t_gIndex,
        std::string& t_ref,
        std::string& t_alt,
        std::string& t_marker,
        uint32_t& t_pd,
        std::string& t_chr,
        double& t_altFreq,
        double& t_altCounts,
        double& t_missingRate,
        double& t_imputeInfo,
        bool& t_isOutputIndexForMissing,
        std::vector<uint>& t_indexForMissing,
        bool& t_isOnlyOutputNonZero,
        std::vector<uint>& t_indexForNonZero,
        bool& t_isBoolRead,
        arma::vec& t_GVec,
        arma::mat& t_GByAncestry,
        bool t_isImputation
    );

    void getOneMarkerAdmixedField(
        uint64_t& t_gIndex_prev,
        uint64_t& t_gIndex,
        std::string& t_ref,
        std::string& t_alt,
        std::string& t_marker,
        uint32_t& t_pd,
        std::string& t_chr,
        double& t_altFreq,
        double& t_altCounts,
        double& t_missingRate,
        double& t_imputeInfo,
        bool& t_isOutputIndexForMissing,
        std::vector<uint>& t_indexForMissing,
        bool& t_isOnlyOutputNonZero,
        std::vector<uint>& t_indexForNonZero,
        bool& t_isBoolRead,
        arma::vec& t_GVec,
        bool t_isImputation,
        const std::string& t_vcfField
    );

    uint32_t getN0() const;
    uint32_t getN() const;
    uint32_t getM0() const;
    uint32_t getM() const;

    void closegenofile();
    void set_iterator(const std::string& chr, int64_t beg_pos, int64_t end_pos);
    void move_forward_iterator(int n_marker);
    bool check_iterator_end();
    uint32_t blockIdForMarker(uint64_t t_gIndex);
    double get_io_seconds() const;
    double get_decode_seconds() const;

private:
    struct AncBlockRecord {
        uint32_t block_id = 0;
        std::string chr;
        int64_t start_pos = 0;
        int64_t end_pos = 0;
        uint64_t anc_offset = 0;
    };

    std::string prefix_;
    Meta* meta_ = nullptr;
    MksStream* common_stream_ = nullptr;
    MksStream* rare_stream_ = nullptr;

    FILE* common_bin_ = nullptr;
    FILE* rare_bin_ = nullptr;
    FILE* anc_bin_ = nullptr;
    FILE* anc_idx_ = nullptr;

    std::vector<std::string> samples_;
    std::vector<int32_t> model_to_raw_sample_;
    std::vector<int32_t> raw_sample_to_model_;

    std::vector<uint64_t> common_words_;
    std::vector<uint8_t> hap_ancestry_;
    std::vector<AncBlockRecord> anc_blocks_;
    std::vector<RareCarrierPacked> cached_rare_carriers_;
    arma::vec cached_total_dosage_;
    arma::mat cached_dosage_by_ancestry_;
    arma::mat cached_ancestry_counts_;
    std::vector<uint> cached_total_nonzero_;
    std::vector<std::vector<uint>> cached_dosage_nonzero_by_ancestry_;
    std::vector<std::vector<uint>> cached_ancestry_nonzero_by_ancestry_;
    std::vector<double> cached_dosage_alt_counts_by_ancestry_;
    std::vector<double> cached_ancestry_counts_by_ancestry_;
    double cached_total_alt_counts_ = 0.0;
    uint32_t cached_common_global_variant_index_ = UINT32_MAX;
    uint32_t cached_rare_global_variant_index_ = UINT32_MAX;
    uint32_t cached_dosage_global_variant_index_ = UINT32_MAX;
    uint32_t cached_anc_block_id_ = UINT32_MAX;
    uint32_t cached_ancestry_counts_block_id_ = UINT32_MAX;
    size_t anc_block_cursor_ = 0;
    bool streaming_mode_ = false;
    bool has_region_filter_ = false;
    std::string region_chr_;
    int64_t region_start_ = 1;
    int64_t region_end_ = 0;

    uint64_t n_common_ = 0;
    uint64_t n_rare_ = 0;
    uint64_t common_bin_pos_ = 0;
    uint64_t rare_bin_pos_ = 0;
    uint64_t anc_bin_pos_ = 0;
    uint64_t anc_idx_pos_ = 8;
    double io_seconds_ = 0.0;
    double decode_seconds_ = 0.0;

    const VariantRecord* current_record() const;
    void advance_current();
    bool advance_to(uint32_t global_variant_index);
    void advance_to_next_in_filter();
    const VariantRecord* record_for_request(uint64_t t_gIndex, bool& t_isBoolRead);
    void fill_marker_fields(
        const VariantRecord& record,
        std::string& t_ref,
        std::string& t_alt,
        std::string& t_marker,
        uint32_t& t_pd,
        std::string& t_chr,
        double& t_missingRate,
        double& t_imputeInfo,
        bool t_isImputation
    ) const;
    void read_dosage(
        const VariantRecord& record,
        arma::vec& dosages,
        std::vector<uint>& indexForNonZero,
        double& altCounts
    );
    void ensure_marker_dosage_cache(const VariantRecord& record);
    void ensure_ancestry_count_cache(uint32_t block_id);
    void load_common_words(const VariantRecord& record);
    void load_rare_carriers(const VariantRecord& record);
    void read_dosage_for_ancestry(
        const VariantRecord& record,
        uint32_t ancestry,
        arma::vec& dosages,
        std::vector<uint>& indexForNonZero,
        double& altCounts
    );
    void read_ancestry_count(
        const VariantRecord& record,
        uint32_t ancestry,
        arma::vec& ancestryCounts,
        std::vector<uint>& indexForNonZero,
        double& altCounts
    );
    void read_dosage_by_ancestry(
        const VariantRecord& record,
        arma::vec& dosages,
        arma::mat& dosageByAncestry,
        std::vector<uint>& indexForNonZero,
        double& altCounts
    );
    void read_anc_blocks();
    uint32_t block_id_for_record(const VariantRecord& record);
    void load_ancestry_block(uint32_t block_id);
    void tracked_seek(FILE* fp, uint64_t offset, const std::string& path, uint64_t& current_offset);
    void tracked_read_exact(FILE* fp, void* data, size_t bytes, const std::string& path, uint64_t& current_offset);
    size_t tracked_fread(FILE* fp, void* data, size_t size, size_t count, const std::string& path, uint64_t& current_offset);
    uint32_t tracked_read_u32_le(FILE* fp, const std::string& path, uint64_t& current_offset);
    uint64_t tracked_read_u64_le(FILE* fp, const std::string& path, uint64_t& current_offset);
};

} // namespace TRACTOR_HYBRID

#endif
