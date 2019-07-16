#!/bin/bash

# Define directory structure; for scripts
ANA_SRC=scripts/analysis
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

# Define directory structure; for output
MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

# Define directory structure; for intermediates
DATA=data
LOG_ANA=log/analysis
LOG_FIG=log/figures
LOG_TAB=log/tables

# Define commands
RDAT=R CMD BATCH --vanilla $< $(LOG_ANA)/$(<F).Rout
RFIG=R CMD BATCH --vanilla $< $(LOG_FIG)/$(<F).Rout
RTAB=R CMD BATCH --vanilla $< $(LOG_TAB)/$(<F).Rout

# All
all: ## Run all parts of the makefile
all: analysis figures tables clean

# Directories
dir_data: ## Make data directory
dir_data:
	test ! -d $(DATA) && mkdir $(DATA) || exit 0
dir_manuscript: ## Make manuscript directory tree
dir_manuscript:
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || exit 0
	test ! -d $(TAB_DIR) && mkdir $(TAB_DIR) || exit 0
	test ! -d $(FIG_DIR) && mkdir $(FIG_DIR) || exit 0
dir_logs: ## Make logs directory tree
dir_logs:
	test ! -d $(LOG_ANA) && mkdir $(LOG_ANA) || exit 0
	test ! -d $(LOG_FIG) && mkdir $(LOG_FIG) || exit 0
	test ! -d $(LOG_TAB) && mkdir $(LOG_TAB) || exit 0
	
analysis: ## Run the analysis
analysis: dir_data \
	dir_logs \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds \
	$(DATA)/gene_counts.rds	\
	$(DATA)/transformed_counts.rds \
	$(DATA)/peak_counts.rds \
	$(DATA)/binding_data.rds \
	$(DATA)/factor_occupancy.rds \
	$(DATA)/histone_modification.rds \
	$(DATA)/hm_occupancy.rds \
	$(DATA)/ddcor.rds \
	$(DATA)/dgca.rds \
	$(DATA)/interactions.rds \
	$(DATA)/factor_targets.rds \
	$(DATA)/deg_res.rds \
	$(DATA)/dep_res.rds \
	$(DATA)/occupancy_res.rds \
	$(DATA)/factor_track_signal.rds \
	$(DATA)/hm_track_signal.rds \
	$(DATA)/polr2a_track_signal.rds \
	$(DATA)/data_tracks.rds \
	$(DATA)/data_tracks_tidy.rds \
	$(DATA)/data_tracks_tissue.rds \
	$(DATA)/data_tracks_tissue_tidy.rds \
	$(DATA)/coverage_tracks.rds \
	$(DATA)/coverage_tracks_diff.rds \
	$(DATA)/peak_overlaps.rds

figures: ## Generate the figures
figures: dir_manuscript \
	dir_logs \
	$(FIG_DIR)/markers.png \
	$(FIG_DIR)/mds_gene_group.png \
	$(FIG_DIR)/mds_binding_factors.png \
	$(FIG_DIR)/volcanos.png \
	$(FIG_DIR)/volcanos_binding.png \
	$(FIG_DIR)/correlations_factor_stage.png \
	$(FIG_DIR)/occupancy_factor_direction.png \
	$(FIG_DIR)/modification_factor_direction.png \
	$(FIG_DIR)/transcription_factor_direction.png \
	$(FIG_DIR)/string_chip_networks.png \
	$(FIG_DIR)/profiles_autophagy_genes.png \
	$(FIG_DIR)/profiles_autophagy_tfs.png \
	$(FIG_DIR)/profiles_adipogenic_tfs.png \
	$(FIG_DIR)/peaks_autophagy_genes.png \
	$(FIG_DIR)/peaks_autophagy_tfs.png \
	$(FIG_DIR)/peaks_adipogenic_tfs.png \
	$(FIG_DIR)/coexpres_adipogenic_autophagy_tf.png \
	$(FIG_DIR)/coexpres_adipogenic_autophagy_genes.png \
	$(FIG_DIR)/coexpres_adipogenic_adipogenic.png \
	$(FIG_DIR)/factor_correlations_heatmap.png \
	$(FIG_DIR)/annotation_correlations_heatmap.png \
	$(FIG_DIR)/factor_cofac_hm.png \
	$(FIG_DIR)/tf_tracks.png \
	$(FIG_DIR)/hm_tracks.png \
	$(FIG_DIR)/pparg_signal.png \
	$(FIG_DIR)/cebpb_signal.png \
	$(FIG_DIR)/med1_signal.png \
	$(FIG_DIR)/rxrg_signal.png \
	$(FIG_DIR)/cebpb_signal_direct.png \
	$(FIG_DIR)/pparg_signal_direct.png \
	$(FIG_DIR)/pparg_signal_tissue.png \
	$(FIG_DIR)/cebpb_signal_tissue.png \
	$(FIG_DIR)/peak_overlap_heatmap.png \
	$(FIG_DIR)/transcription_expression.png

tables: ## Generate the tables
tables: dir_manuscript \
	dir_logs \
	$(TAB_DIR)/datasets.tex \
	$(TAB_DIR)/deg_markers.tex \
	$(TAB_DIR)/deg_fc.tex \
	$(TAB_DIR)/dep_fc.tex \
	$(TAB_DIR)/dep_fc_genes.tex \
	$(TAB_DIR)/dep_fc_tf.tex \
	$(TAB_DIR)/variance_explained.tex \
	$(TAB_DIR)/variance_explained_chip.tex
	
# Analysis
$(DATA)/go_annotation.rds: $(ANA_SRC)/go_annotation.R
	$(RDAT)
$(DATA)/tf_annotation.rds: $(ANA_SRC)/tf_annotation.R
	$(RDAT)
$(DATA)/gene_counts.rds: $(ANA_SRC)/gene_counts.R
	$(RDAT)
$(DATA)/transformed_counts.rds: $(ANA_SRC)/transformed_counts.R \
	$(DATA)/gene_counts.rds
	$(RDAT)
$(DATA)/peak_counts.rds: $(ANA_SRC)/peak_counts.R \
	$(DATA)/peak_counts.rda
	$(RDAT)
$(DATA)/binding_data.rds: $(ANA_SRC)/binding_data.R \
	$(DATA)/peak_counts.rds
	$(RDAT)
$(DATA)/factor_occupancy.rds: $(ANA_SRC)/factor_occupancy.R \
	$(DATA)/binding_data.rds
	$(RDAT)
$(DATA)/histone_modification.rds: $(ANA_SRC)/histone_modification.R \
	$(DATA)/peak_counts.rds
	$(RDAT)
$(DATA)/hm_occupancy.rds: $(ANA_SRC)/hm_occupancy.R \
	$(DATA)/histone_modification.rds
	$(RDAT)
$(DATA)/ddcor.rds: $(ANA_SRC)/ddcor.R \
	$(DATA)/transformed_counts.rds \
	$(DATA)/go_annotation.rds
	$(RDAT)
$(DATA)/dgca.rds: $(ANA_SRC)/dgca.R \
	$(DATA)/ddcor.rds
	$(RDAT)
$(DATA)/interactions.rds: $(ANA_SRC)/interactions.R \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RDAT)
$(DATA)/factor_targets.rds: $(ANA_SRC)/factor_targets.R \
	$(DATA)/peak_counts.rds
	$(RDAT)
$(DATA)/deg_res.rds: $(ANA_SRC)/deg_res.R \
	$(DATA)/gene_counts.rds
	$(RDAT)
$(DATA)/dep_res.rds: $(ANA_SRC)/dep_res.R \
	$(DATA)/binding_data.rds \
	$(DATA)/histone_modification.rds
	$(RDAT)
$(DATA)/occupancy_res.rds: $(ANA_SRC)/occupancy_res.R \
	$(DATA)/factor_occupancy.rds
	$(RDAT)
$(DATA)/factor_track_signal.rds: $(ANA_SRC)/factor_track_signal.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/go_annotation.rds
	$(RDAT)
$(DATA)/hm_track_signal.rds: $(ANA_SRC)/hm_track_signal.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/go_annotation.rds
	$(RDAT)
$(DATA)/polr2a_track_signal.rds: $(ANA_SRC)/polr2a_track_signal.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/go_annotation.rds
	$(RDAT)
$(DATA)/data_tracks.rds: $(ANA_SRC)/data_tracks.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/factor_track_signal.rds
	$(RDAT)
$(DATA)/data_tracks_tidy.rds: $(ANA_SRC)/data_tracks_tidy.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/factor_track_signal.rds \
	$(DATA)/go_annotation.rds
	$(RDAT)
$(DATA)/coverage_tracks.rds: $(ANA_SRC)/coverage_tracks.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/factor_track_signal.rds \
	$(DATA)/hm_track_signal.rds \
	$(DATA)/go_annotation.rds
	$(RDAT)
$(DATA)/coverage_tracks_diff.rds: $(ANA_SRC)/coverage_tracks_diff.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/factor_track_signal.rds \
	$(DATA)/hm_track_signal.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/deg_res.rds
	$(RDAT)
$(DATA)/data_tracks_tissue.rds: $(ANA_SRC)/data_tracks_tissue.R \
	$(DATA)/go_annotation.rds
	$(RDAT)
$(DATA)/data_tracks_tissue_tidy.rds: $(ANA_SRC)/data_tracks_tissue_tidy.R \
	$(DATA)/go_annotation.rds
	$(RDAT)
$(DATA)/peak_overlaps.rds: $(ANA_SRC)/peak_overlaps.R \
	$(DATA)/peak_counts.rds
	$(RDAT)

# Figures
$(FIG_DIR)/markers.png: $(FIG_SRC)/markers.R \
	$(DATA)/gene_counts.rds
	$(RFIG)
$(FIG_DIR)/mds_gene_group.png: $(FIG_SRC)/mds_gene_group.R \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds \
	$(DATA)/transformed_counts.rds
	$(RFIG)
$(FIG_DIR)/mds_binding_factors.png: $(FIG_SRC)/mds_binding_factors.R \
	$(DATA)/binding_data.rds
	$(RFIG)
$(FIG_DIR)/volcanos.png: $(FIG_SRC)/volcanos.R \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds \
	$(DATA)/deg_res.rds
	$(RFIG)
$(FIG_DIR)/volcanos_binding.png: $(FIG_SRC)/volcanos_binding.R \
	$(DATA)/go_annotation.rds \
	$(DATA)/binding_data.rds \
	$(DATA)/dep_res.rds
	$(RFIG)
$(FIG_DIR)/correlations_factor_stage.png: $(FIG_SRC)/correlations_factor_stage.R \
	$(DATA)/deg_res.rds \
	$(DATA)/ddcor.rds
	$(RFIG)
$(FIG_DIR)/occupancy_factor_direction.png: $(FIG_SRC)/occupancy_factor_direction.R \
	$(DATA)/deg_res.rds \
	$(DATA)/factor_occupancy.rds
	$(RFIG)	
$(FIG_DIR)/modification_factor_direction.png: $(FIG_SRC)/modification_factor_direction.R \
	$(DATA)/deg_res.rds \
	$(DATA)/hm_occupancy.rds
	$(RFIG)
$(FIG_DIR)/transcription_factor_direction.png: $(FIG_SRC)/transcription_factor_direction.R \
	$(DATA)/deg_res.rds \
	$(DATA)/factor_occupancy.rds
	$(RFIG)	
$(FIG_DIR)/string_chip_networks.png: $(FIG_SRC)/string_chip_networks.R \
	$(DATA)/interactions.rds \
	$(DATA)/factor_targets.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RFIG)	
$(FIG_DIR)/profiles_autophagy_genes.png: $(FIG_SRC)/profiles_autophagy_genes.R \
	$(DATA)/gene_counts.rds \
	$(DATA)/factor_occupancy.rds
	$(RFIG)
$(FIG_DIR)/profiles_autophagy_tfs.png: $(FIG_SRC)/profiles_autophagy_tfs.R \
	$(DATA)/gene_counts.rds \
	$(DATA)/factor_occupancy.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RFIG)
$(FIG_DIR)/profiles_adipogenic_tfs.png: $(FIG_SRC)/profiles_adipogenic_tfs.R \
	$(DATA)/gene_counts.rds \
	$(DATA)/factor_occupancy.rds
	$(RFIG)
$(FIG_DIR)/peaks_autophagy_genes.png: $(FIG_SRC)/peaks_autophagy_genes.R \
	$(DATA)/binding_data.rds
	$(RFIG)
$(FIG_DIR)/peaks_autophagy_tfs.png: $(FIG_SRC)/peaks_autophagy_tfs.R \
	$(DATA)/binding_data.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RFIG)
$(FIG_DIR)/peaks_adipogenic_tfs.png: $(FIG_SRC)/peaks_adipogenic_tfs.R \
	$(DATA)/binding_data.rds
	$(RFIG)
$(FIG_DIR)/coexpres_adipogenic_autophagy_tf.png: $(FIG_SRC)/coexpres_adipogenic_autophagy_tf.R \
	$(DATA)/dgca.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RFIG)
$(FIG_DIR)/coexpres_adipogenic_autophagy_genes.png: $(FIG_SRC)/coexpres_adipogenic_autophagy_genes.R \
	$(DATA)/dgca.rds
	$(RFIG)
$(FIG_DIR)/coexpres_adipogenic_adipogenic.png: $(FIG_SRC)/coexpres_adipogenic_adipogenic.R \
	$(DATA)/dgca.rds
	$(RFIG)
$(FIG_DIR)/factor_correlations_heatmap.png: $(FIG_SRC)/factor_correlations_heatmap.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/go_annotation.rds
	$(RFIG)
$(FIG_DIR)/annotation_correlations_heatmap.png: $(FIG_SRC)/annotation_correlations_heatmap.R \
	$(DATA)/peak_counts.rds $(DATA)/go_annotation.rds
	$(RFIG)
$(FIG_DIR)/factor_cofac_hm.png: $(FIG_SRC)/factor_cofac_hm.R \
	$(DATA)/peak_counts.rds $(DATA)/go_annotation.rds
	$(RFIG)
$(FIG_DIR)/tf_tracks.png: $(FIG_SRC)/tf_tracks.R \
	$(DATA)/coverage_tracks.rds
	$(RFIG)
$(FIG_DIR)/hm_tracks.png: $(FIG_SRC)/hm_tracks.R \
	$(DATA)/coverage_tracks.rds
	$(RFIG)
$(FIG_DIR)/pparg_signal.png: $(FIG_SRC)/pparg_signal.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/cebpb_signal.png: $(FIG_SRC)/cebpb_signal.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/med1_signal.png: $(FIG_SRC)/med1_signal.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/rxrg_signal.png: $(FIG_SRC)/rxrg_signal.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/pparg_signal_direct.png: $(FIG_SRC)/pparg_signal_direct.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/cebpb_signal_direct.png: $(FIG_SRC)/cebpb_signal_direct.R \
	$(DATA)/data_tracks.rds
	$(RFIG)
$(FIG_DIR)/cebpb_signal_tissue.png: $(FIG_SRC)/cebpb_signal_tissue.R \
	$(DATA)/data_tracks_tissue.rds
	$(RFIG)
$(FIG_DIR)/pparg_signal_tissue.png: $(FIG_SRC)/pparg_signal_tissue.R \
	$(DATA)/data_tracks_tissue.rds
	$(RFIG)
$(FIG_DIR)/transcription_expression.png: $(FIG_SRC)/transcription_expression.R \
	$(DATA)/transformed_counts.rds $(DATA)/go_annotation.rds \
	$(DATA)/factor_occupancy.rds
	$(RFIG)
$(FIG_DIR)/peak_overlap_heatmap.png: $(FIG_SRC)/peak_overlap_heatmap.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/peak_overlaps.rds
	$(RFIG)

# Tables
$(TAB_DIR)/datasets.tex: $(TAB_SRC)/datasets.R \
	$(DATA)/peak_counts.rds \
	$(DATA)/gene_counts.rds
	$(RTAB)
$(TAB_DIR)/deg_markers.tex: $(TAB_SRC)/deg_markers.R \
	$(DATA)/deg_res.rds
	$(RTAB)
$(TAB_DIR)/deg_fc.tex: $(TAB_SRC)/deg_fc.R \
	$(DATA)/deg_res.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RTAB)
$(TAB_DIR)/dep_fc.tex: $(TAB_SRC)/dep_fc.R \
	$(DATA)/binding_data.rds \
	$(DATA)/dep_res.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RTAB)
$(TAB_DIR)/dep_fc_genes.tex: $(TAB_SRC)/dep_fc_genes.R \
	$(DATA)/binding_data.rds \
	$(DATA)/dep_res.rds
	$(RTAB)
$(TAB_DIR)/dep_fc_tf.tex: $(TAB_SRC)/dep_fc_tf.R \
	$(DATA)/binding_data.rds \
	$(DATA)/dep_res.rds
	$(RTAB)
$(TAB_DIR)/variance_explained.tex: $(TAB_SRC)/variance_explained.R \
	$(DATA)/transformed_counts.rds \
	$(DATA)/go_annotation.rds \
	$(DATA)/tf_annotation.rds
	$(RTAB)
$(TAB_DIR)/variance_explained_chip.tex: $(TAB_SRC)/variance_explained_chip.R $(DATA)/binding_data.rds
	$(RTAB)

# Clean Up
.PHONY: clean
clean: ## Clean up
clean:
	rm -f *.pdf
	rm -f *.RData

# Source: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help
help: ## Print the current page
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help
