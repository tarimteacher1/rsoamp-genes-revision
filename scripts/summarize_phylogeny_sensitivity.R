#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) >= 1) args[[1]] else "/path/to/rsoamp/revision_R1_20260902"
phylo_dir <- file.path(base, "02_phylogeny")

families <- c("Defensin", "Snakin_GASA", "nsLTP")
regions <- c("full", "core")
trims <- c("untrimmed", "gappyout")

parse_support <- function(label) {
  if (!length(label) || is.na(label) || !nzchar(label)) return(c(NA_real_, NA_real_))
  parts <- strsplit(label, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(c(NA_real_, NA_real_))
  suppressWarnings(as.numeric(parts))
}

descendant_tips <- function(tree, node) {
  ntip <- length(tree$tip.label)
  children <- tree$edge[tree$edge[, 1] == node, 2]
  tips <- children[children <= ntip]
  internal <- children[children > ntip]
  if (length(internal)) {
    tips <- c(tips, unlist(lapply(internal, function(x) descendant_tips(tree, x))))
  }
  unique(tips)
}

canonical_side <- function(a, b) {
  if (length(a) < length(b)) return(sort(a))
  if (length(b) < length(a)) return(sort(b))
  sa <- paste(sort(a), collapse = ";")
  sb <- paste(sort(b), collapse = ";")
  if (sa <= sb) sort(a) else sort(b)
}

branch_rows <- list()
row_index <- 1L

for (family in families) {
  for (region in regions) {
    for (trim in trims) {
      strategy <- paste(region, trim, sep = ".")
      tree_path <- file.path(phylo_dir, family, sprintf("%s.%s.%s.treefile", family, region, trim))
      if (!file.exists(tree_path) || file.info(tree_path)$size == 0) {
        stop(sprintf("Missing tree: %s", tree_path))
      }
      tree <- read.tree(tree_path)
      ntip <- length(tree$tip.label)
      all_tips <- sort(tree$tip.label)

      for (child in unique(tree$edge[, 2])) {
        if (child <= ntip) next
        side_a <- sort(tree$tip.label[descendant_tips(tree, child)])
        side_b <- setdiff(all_tips, side_a)
        split_side <- canonical_side(side_a, side_b)
        node_label <- tree$node.label[[child - ntip]]
        support <- parse_support(node_label)

        rso_a <- sub("^RSO_", "", side_a[grepl("^RSO_", side_a)])
        rso_b <- sub("^RSO_", "", side_b[grepl("^RSO_", side_b)])
        rso_group <- if (!length(rso_a) || !length(rso_b)) {
          character()
        } else {
          canonical_side(rso_a, rso_b)
        }

        branch_rows[[row_index]] <- data.frame(
          family = family,
          strategy = strategy,
          alignment_region = region,
          trimming = trim,
          node = child,
          sh_alrt = support[[1]],
          ufboot = support[[2]],
          joint_support = !is.na(support[[1]]) && !is.na(support[[2]]) && support[[1]] >= 80 && support[[2]] >= 95,
          side_a_taxa = paste(side_a, collapse = ";"),
          side_b_taxa = paste(side_b, collapse = ";"),
          canonical_unrooted_split = paste(split_side, collapse = ";"),
          canonical_split_size = length(split_side),
          rso_group = paste(sort(rso_group), collapse = ";"),
          rso_group_size = length(rso_group),
          tree_file = tree_path,
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }
  }
}

branches <- do.call(rbind, branch_rows)
write.table(
  branches,
  file.path(phylo_dir, "phylogeny_branch_support.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

split_taxa <- function(text) {
  if (is.na(text) || !nzchar(text)) character() else strsplit(text, ";", fixed = TRUE)[[1]]
}

reference_class <- function(taxa) {
  counts <- c(
    positive_nsLTP = sum(grepl("^NSPOS_", taxa)),
    explicit_2S = sum(grepl("^NEG2S_", taxa)),
    ambiguous_prolamin = sum(grepl("^AMBPRO_", taxa)),
    tamarix = sum(grepl("^TAU_", taxa))
  )
  active <- names(counts)[counts > 0]
  call <- if (!length(active)) {
    "NO_LABELED_REFERENCE"
  } else if (identical(active, "positive_nsLTP")) {
    "POSITIVE_NSLTP_ONLY"
  } else if ("explicit_2S" %in% active) {
    "EXPLICIT_2S_PRESENT"
  } else if ("ambiguous_prolamin" %in% active) {
    "AMBIGUOUS_PROLAMIN_PRESENT"
  } else if (identical(active, "tamarix")) {
    "TAMARIX_ONLY"
  } else {
    "MIXED_REFERENCE_CLASSES"
  }
  c(counts, call = call)
}

ns_branches <- branches[branches$family == "nsLTP", ]
ns_taxa <- unique(unlist(lapply(c(ns_branches$side_a_taxa, ns_branches$side_b_taxa), split_taxa)))
ns_rso_taxa <- sort(ns_taxa[grepl("^RSO_", ns_taxa)])
neighborhood_rows <- list()
neighborhood_index <- 1L
for (strategy in sort(unique(ns_branches$strategy))) {
  strategy_edges <- ns_branches[ns_branches$strategy == strategy & ns_branches$joint_support, ]
  for (taxon in ns_rso_taxa) {
    candidates <- list()
    candidate_index <- 1L
    for (edge_index in seq_len(nrow(strategy_edges))) {
      edge <- strategy_edges[edge_index, ]
      side_a <- split_taxa(edge$side_a_taxa)
      side_b <- split_taxa(edge$side_b_taxa)
      side <- if (taxon %in% side_a) side_a else side_b
      classes <- reference_class(side)
      if (classes[["call"]] == "NO_LABELED_REFERENCE") next
      candidates[[candidate_index]] <- data.frame(
        family = "nsLTP",
        strategy = strategy,
        candidate = sub("^RSO_", "", taxon),
        supported_neighborhood_size = length(side),
        rso_members = paste(sort(sub("^RSO_", "", side[grepl("^RSO_", side)])), collapse = ";"),
        positive_nsLTP_references = as.integer(classes[["positive_nsLTP"]]),
        explicit_2S_references = as.integer(classes[["explicit_2S"]]),
        ambiguous_prolamin_references = as.integer(classes[["ambiguous_prolamin"]]),
        tamarix_references = as.integer(classes[["tamarix"]]),
        reference_neighborhood_call = classes[["call"]],
        sh_alrt = edge$sh_alrt,
        ufboot = edge$ufboot,
        stringsAsFactors = FALSE
      )
      candidate_index <- candidate_index + 1L
    }
    if (length(candidates)) {
      candidate_table <- do.call(rbind, candidates)
      candidate_table <- candidate_table[order(candidate_table$supported_neighborhood_size, -candidate_table$sh_alrt, -candidate_table$ufboot), ]
      neighborhood_rows[[neighborhood_index]] <- candidate_table[1, ]
    } else {
      neighborhood_rows[[neighborhood_index]] <- data.frame(
        family = "nsLTP", strategy = strategy, candidate = sub("^RSO_", "", taxon),
        supported_neighborhood_size = NA_integer_, rso_members = "",
        positive_nsLTP_references = 0L, explicit_2S_references = 0L,
        ambiguous_prolamin_references = 0L, tamarix_references = 0L,
        reference_neighborhood_call = "NO_SUPPORTED_REFERENCE_NEIGHBORHOOD",
        sh_alrt = NA_real_, ufboot = NA_real_, stringsAsFactors = FALSE
      )
    }
    neighborhood_index <- neighborhood_index + 1L
  }
}
neighborhood_table <- do.call(rbind, neighborhood_rows)
write.table(
  neighborhood_table,
  file.path(phylo_dir, "nsltp_supported_reference_neighborhoods.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

candidate_neighborhoods <- split(neighborhood_table, neighborhood_table$candidate)
neighborhood_summary <- do.call(rbind, lapply(candidate_neighborhoods, function(x) {
  calls <- table(factor(
    x$reference_neighborhood_call,
    levels = c("POSITIVE_NSLTP_ONLY", "EXPLICIT_2S_PRESENT", "AMBIGUOUS_PROLAMIN_PRESENT", "TAMARIX_ONLY", "MIXED_REFERENCE_CLASSES", "NO_SUPPORTED_REFERENCE_NEIGHBORHOOD")
  ))
  data.frame(
    candidate = x$candidate[[1]],
    strategies_positive_nsLTP_only = calls[["POSITIVE_NSLTP_ONLY"]],
    strategies_explicit_2S_present = calls[["EXPLICIT_2S_PRESENT"]],
    strategies_ambiguous_prolamin_present = calls[["AMBIGUOUS_PROLAMIN_PRESENT"]],
    strategies_tamarix_only = calls[["TAMARIX_ONLY"]],
    strategies_mixed_reference_classes = calls[["MIXED_REFERENCE_CLASSES"]],
    strategies_without_supported_reference_neighborhood = calls[["NO_SUPPORTED_REFERENCE_NEIGHBORHOOD"]],
    phylogenetic_reference_interpretation = if (calls[["POSITIVE_NSLTP_ONLY"]] >= 3) {
      "CONSISTENT_POSITIVE_NSLTP_NEIGHBORHOOD"
    } else if (calls[["EXPLICIT_2S_PRESENT"]] >= 2) {
      "REPEATED_EXPLICIT_2S_ASSOCIATION"
    } else if (calls[["AMBIGUOUS_PROLAMIN_PRESENT"]] >= 2 || calls[["MIXED_REFERENCE_CLASSES"]] >= 2) {
      "REPEATED_AMBIGUOUS_OR_MIXED_ASSOCIATION"
    } else {
      "PHYLOGENETIC_NEIGHBORHOOD_UNRESOLVED"
    },
    stringsAsFactors = FALSE
  )
}))
write.table(
  neighborhood_summary,
  file.path(phylo_dir, "nsltp_reference_neighborhood_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

supported <- branches[branches$joint_support & branches$rso_group_size >= 1 & nzchar(branches$rso_group), ]
if (nrow(supported)) {
  key <- paste(supported$family, supported$rso_group, sep = "||")
  groups <- split(supported, key)
  group_rows <- lapply(groups, function(x) {
    best_by_strategy <- aggregate(
      cbind(sh_alrt, ufboot) ~ strategy,
      data = x,
      FUN = max,
      na.rm = TRUE
    )
    n_strategies <- length(unique(x$strategy))
    data.frame(
      family = x$family[[1]],
      rso_group = x$rso_group[[1]],
      rso_group_size = x$rso_group_size[[1]],
      supported_strategy_count = n_strategies,
      supported_strategies = paste(sort(unique(x$strategy)), collapse = ";"),
      support_values = paste(
        sprintf("%s=%.1f/%.1f", best_by_strategy$strategy, best_by_strategy$sh_alrt, best_by_strategy$ufboot),
        collapse = ";"
      ),
      stability_class = if (n_strategies == 4) {
        "STABLE_4_OF_4"
      } else if (n_strategies == 3) {
        "STABLE_3_OF_4"
      } else if (n_strategies == 2) {
        "SENSITIVE_2_OF_4"
      } else {
        "SINGLE_STRATEGY_ONLY"
      },
      stringsAsFactors = FALSE
    )
  })
  group_table <- do.call(rbind, group_rows)
  group_table <- group_table[order(group_table$family, -group_table$supported_strategy_count, group_table$rso_group_size, group_table$rso_group), ]
} else {
  group_table <- data.frame(
    family = character(), rso_group = character(), rso_group_size = integer(),
    supported_strategy_count = integer(), supported_strategies = character(),
    support_values = character(), stability_class = character()
  )
}

write.table(
  group_table,
  file.path(phylo_dir, "rso_supported_group_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

all_rso <- do.call(rbind, lapply(families, function(family) {
  manifest <- read.delim(file.path(phylo_dir, family, sprintf("%s.taxon_manifest.tsv", family)), stringsAsFactors = FALSE)
  manifest <- manifest[grepl("^RSO_", manifest$taxon_id), , drop = FALSE]
  data.frame(family = family, member_id = sub("^RSO_", "", manifest$taxon_id), stringsAsFactors = FALSE)
}))
all_rso <- unique(all_rso)
taxon_rows <- lapply(seq_len(nrow(all_rso)), function(index) {
  taxon <- all_rso$member_id[[index]]
  family <- all_rso$family[[index]]
  family_groups <- group_table[group_table$family == family, , drop = FALSE]
  in_group <- if (nrow(family_groups)) {
    vapply(strsplit(family_groups$rso_group, ";", fixed = TRUE), function(x) taxon %in% x, logical(1))
  } else {
    logical()
  }
  stable <- in_group & family_groups$supported_strategy_count >= 3 & family_groups$rso_group_size >= 2
  data.frame(
    member_id = taxon,
    family = family,
    stable_supported_group_count = sum(stable),
    smallest_stable_group_size = if (any(stable)) min(family_groups$rso_group_size[stable]) else NA_integer_,
    placement_status = if (any(stable)) "STABLE_SUPPORTED_GROUP" else "UNRESOLVED_ACROSS_SENSITIVITY_ANALYSES",
    stringsAsFactors = FALSE
  )
})
taxon_table <- if (length(taxon_rows)) do.call(rbind, taxon_rows) else data.frame()
write.table(
  taxon_table,
  file.path(phylo_dir, "rso_phylogenetic_placement_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

summary_rows <- do.call(rbind, lapply(families, function(family) {
  subset <- branches[branches$family == family, ]
  data.frame(
    family = family,
    total_internal_edges_across_four_strategies = nrow(subset),
    jointly_supported_edges_across_four_strategies = sum(subset$joint_support),
    stable_rso_groups_4_of_4 = sum(group_table$family == family & group_table$supported_strategy_count == 4 & group_table$rso_group_size >= 2),
    stable_rso_groups_3_of_4 = sum(group_table$family == family & group_table$supported_strategy_count == 3 & group_table$rso_group_size >= 2),
    interpretation = "Supported clades are sensitivity-qualified; unsupported backbone placement is not assigned to a formal subfamily",
    stringsAsFactors = FALSE
  )
}))
write.table(
  summary_rows,
  file.path(phylo_dir, "phylogeny_sensitivity_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

alignment_rows <- list()
alignment_index <- 1L
for (family in families) {
  for (region in regions) {
    for (trim in trims) {
      strategy <- paste(region, trim, sep = ".")
      report_path <- file.path(phylo_dir, family, sprintf("%s.%s.%s.iqtree", family, region, trim))
      report <- readLines(report_path, warn = FALSE)
      input_line <- grep("^Input data:", report, value = TRUE)[1]
      model_line <- grep("^Best-fit model according to BIC:", report, value = TRUE)[1]
      likelihood_line <- grep("^Log-likelihood of the tree:", report, value = TRUE)[1]
      input_match <- regexec("Input data: ([0-9]+) sequences with ([0-9]+) amino-acid sites", input_line)
      input_values <- regmatches(input_line, input_match)[[1]]
      likelihood_match <- regexec("Log-likelihood of the tree: (-?[0-9.]+)", likelihood_line)
      likelihood_values <- regmatches(likelihood_line, likelihood_match)[[1]]
      if (length(input_values) != 3 || length(likelihood_values) != 2) {
        stop(sprintf("Could not parse IQ-TREE report: %s", report_path))
      }
      alignment_rows[[alignment_index]] <- data.frame(
        family = family,
        strategy = strategy,
        alignment_region = region,
        trimming = trim,
        sequence_count = as.integer(input_values[[2]]),
        alignment_sites = as.integer(input_values[[3]]),
        best_fit_model_BIC = sub("^Best-fit model according to BIC: ", "", model_line),
        tree_log_likelihood = as.numeric(likelihood_values[[2]]),
        iqtree_report = report_path,
        tree_file = file.path(phylo_dir, family, sprintf("%s.%s.%s.treefile", family, region, trim)),
        stringsAsFactors = FALSE
      )
      alignment_index <- alignment_index + 1L
    }
  }
}
alignment_table <- do.call(rbind, alignment_rows)
write.table(
  alignment_table,
  file.path(phylo_dir, "alignment_comparison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary_markdown <- c(
  "# Phylogenetic sensitivity summary",
  "",
  "Status: PASS",
  "",
  "Each family was analyzed under four pre-specified strategies: full precursor or mature/core region, each with an untrimmed or trimAl `gappyout` alignment. IQ-TREE model selection used the same candidate model set, and branch support was evaluated with 1,000 SH-aLRT and 1,000 ultrafast bootstrap replicates.",
  "",
  "A branch was jointly supported only when SH-aLRT was at least 80 and UFBoot was at least 95. A sensitivity-stable R. soongarica group required the same multi-member bipartition in at least three of four strategies.",
  "",
  sprintf(
    "The four-strategy analyses identified %d, %d, and %d multi-member groups stable in all four strategies for Defensin, Snakin/GASA, and nsLTP, respectively; nsLTP additionally contained %d group stable in three strategies.",
    summary_rows$stable_rso_groups_4_of_4[summary_rows$family == "Defensin"],
    summary_rows$stable_rso_groups_4_of_4[summary_rows$family == "Snakin_GASA"],
    summary_rows$stable_rso_groups_4_of_4[summary_rows$family == "nsLTP"],
    summary_rows$stable_rso_groups_3_of_4[summary_rows$family == "nsLTP"]
  ),
  "",
  "These groups are reported as sensitivity-qualified supported clades rather than formal subfamilies. Members without a stable multi-member placement remain unresolved. Because the proteins are short and often highly divergent, unresolved backbone relationships are treated as a data limitation rather than manually forced into bifurcating groups.",
  "",
  "Complete, uncollapsed trees and the branch-level support table are retained for reproducibility; low-support nodes may be collapsed only in the display copy to show unresolved polytomies."
)
writeLines(summary_markdown, file.path(phylo_dir, "phylogeny_sensitivity_summary.md"))

writeLines(
  c(
    "PASS",
    sprintf("Parsed %d internal branches across 12 tree analyses", nrow(branches)),
    "Joint support threshold: SH-aLRT >= 80 and UFBoot >= 95",
    "Stable groups require the same R. soongarica bipartition in at least three of four alignment strategies"
  ),
  file.path(phylo_dir, "PHYLOGENY_SENSITIVITY_SUMMARY_COMPLETE.PASS")
)

print(summary_rows)
