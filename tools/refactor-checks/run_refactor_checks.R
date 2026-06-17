#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
update_baseline <- "--update-baseline" %in% args ||
	identical(tolower(Sys.getenv("UPDATE_BASELINE")), "true")
skip_install <- "--no-install" %in% args ||
	identical(tolower(Sys.getenv("DEPMIX_REF_CHECK_INSTALL")), "false")

root <- normalizePath(getwd(), mustWork = TRUE)
desc <- file.path(root, "DESCRIPTION")
if(!file.exists(desc)) {
	stop("Run this script from the depmixS4 package root.", call. = FALSE)
}

script <- file.path(root, "tools", "refactor-checks", "run_refactor_checks.R")
if(!file.exists(script)) {
	stop("Cannot find tools/refactor-checks/run_refactor_checks.R", call. = FALSE)
}

result_root <- file.path(root, "tools", "refactor-checks", "results")
baseline_dir <- file.path(result_root, "baseline")
current_dir <- file.path(result_root, "current")

if(!skip_install && !identical(Sys.getenv("DEPMIX_REF_CHECK_CHILD"), "1")) {
	lib <- file.path(tempdir(), "depmixS4-refactor-lib")
	if(dir.exists(lib)) unlink(lib, recursive = TRUE)
	dir.create(lib, recursive = TRUE, showWarnings = FALSE)
	message("Installing current source into temporary library: ", lib)
	install_cmd <- file.path(R.home("bin"), "R")
	install_status <- system2(install_cmd,
		c("CMD", "INSTALL", "-l", lib, root))
	if(!identical(install_status, 0L)) {
		stop("Temporary R CMD INSTALL failed.", call. = FALSE)
	}
	child_args <- args[args != "--update-baseline"]
	if(update_baseline) child_args <- c(child_args, "--update-baseline")
	child_args <- child_args[child_args != "--no-install"]
	status <- system2(file.path(R.home("bin"), "Rscript"),
		c(script, child_args, "--no-install"),
		env = c(paste0("R_LIBS=", lib), "DEPMIX_REF_CHECK_CHILD=1"))
	quit(status = status)
}

suppressPackageStartupMessages(require(depmixS4))

dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(current_dir, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
run_dir <- file.path(current_dir, paste0("run-", timestamp))
output_dir <- file.path(run_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_file <- file.path(tempdir(), paste0("depmixS4-refactor-plots-", timestamp, ".pdf"))
grDevices::pdf(plot_file)
on.exit({
	while(grDevices::dev.cur() > 1) grDevices::dev.off()
}, add = TRUE)

safe_name <- function(x) {
	gsub("[^A-Za-z0-9_.-]+", "_", x)
}

write_text <- function(path, lines) {
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	writeLines(lines, path, useBytes = TRUE)
}

read_text <- function(path) {
	if(!file.exists(path)) character()
	readLines(path, warn = FALSE)
}

empty_items <- function() {
	data.frame(
		suite = character(),
		item = character(),
		status = character(),
		elapsed_sec = numeric(),
		output_file = character(),
		error = character(),
		stringsAsFactors = FALSE
	)
}

empty_metrics <- function() {
	data.frame(
		suite = character(),
		item = character(),
		metric = character(),
		type = character(),
		value = character(),
		numeric_value = numeric(),
		stringsAsFactors = FALSE
	)
}

metric_rows <- function(suite, item, metrics) {
	if(length(metrics) == 0) return(empty_metrics())
	rows <- lapply(names(metrics), function(name) {
		value <- metrics[[name]]
		if(is.numeric(value) && length(value) == 1 && is.finite(value)) {
			type <- "numeric"
			text <- sprintf("%.15g", value)
			num <- value
		} else {
			type <- class(value)[1]
			if(is.numeric(value)) {
				text <- paste(sprintf("%.15g", value), collapse = ",")
			} else {
				text <- paste(as.character(value), collapse = ",")
			}
			num <- NA_real_
		}
		data.frame(suite = suite, item = item, metric = name, type = type,
			value = text, numeric_value = num, stringsAsFactors = FALSE)
	})
	do.call(rbind, rows)
}

loglik <- function(object) {
	suppressWarnings(unname(as.numeric(logLik(object))))
}

collect_env_metrics <- function(env) {
	metrics <- list()
	for(name in sort(ls(env))) {
		value <- get(name, envir = env)
		if(inherits(value, c("depmix.fitted", "mix.fitted"))) {
			metrics[[paste0("logLik.", name)]] <- loglik(value)
		} else if(is.data.frame(value)) {
			metrics[[paste0("dim.", name)]] <- paste(dim(value), collapse = "x")
			if("state" %in% names(value)) {
				metrics[[paste0("head_state.", name)]] <-
					paste(head(value$state, 10), collapse = ",")
			}
		} else if(is.matrix(value) || is.array(value)) {
			metrics[[paste0("dim.", name)]] <- paste(dim(value), collapse = "x")
		} else if((is.numeric(value) || is.logical(value)) && length(value) <= 10) {
			metrics[[paste0("value.", name)]] <- value
		}
	}
	metrics
}

run_example <- function(topic) {
	env <- new.env(parent = globalenv())
	result <- capture_run(utils::example(topic, package = "depmixS4",
		ask = FALSE, echo = FALSE, local = env, run.dontrun = FALSE,
		character.only = TRUE))
	metrics <- collect_env_metrics(env)
	list(status = if(is.null(result$error)) "OK" else "FAIL",
		elapsed = result$elapsed, output = result$output,
		warnings = result$warnings, error = result$error, metrics = metrics)
}

rd_example_topics <- function() {
	rd_files <- sort(Sys.glob(file.path(root, "man", "*.Rd")))
	topics <- character()
	for(rd in rd_files) {
		lines <- readLines(rd, warn = FALSE)
		if(!any(grepl("^\\\\examples\\{", lines))) next
		name_line <- grep("^\\\\name\\{", lines, value = TRUE)
		if(length(name_line) == 0) next
		topic <- sub("^\\\\name\\{([^}]+)\\}.*$", "\\1", name_line[1])
		topics <- c(topics, topic)
	}
	unique(topics)
}

run_test_file <- function(path) {
	output_path <- tempfile()
	old_wd <- getwd()
	setwd(tempdir())
	on.exit(setwd(old_wd), add = TRUE)
	elapsed <- system.time({
		output <- system2(file.path(R.home("bin"), "Rscript"), normalizePath(path),
			stdout = TRUE, stderr = TRUE)
	})[["elapsed"]]
	status <- attr(output, "status")
	if(is.null(status)) status <- 0L
	write_text(output_path, output)
	list(status = if(identical(status, 0L)) "OK" else "FAIL",
		elapsed = unname(elapsed), output = output, warnings = character(),
		error = if(identical(status, 0L)) "" else paste("exit", status),
		metrics = list(exit_status = status, output_lines = length(output)))
}

capture_run <- function(expr) {
	warnings <- character()
	output <- character()
	error <- NULL
	store <- new.env(parent = emptyenv())
	store$value <- NULL
	elapsed <- system.time({
		output <- capture.output(
			withCallingHandlers(
				tryCatch({
					store$value <- force(expr)
				}, error = function(e) {
					error <<- conditionMessage(e)
					NULL
				}),
				warning = function(w) {
					warnings <<- c(warnings, conditionMessage(w))
					invokeRestart("muffleWarning")
				}
			)
		)
	})[["elapsed"]]
	list(output = output, warnings = warnings, error = error,
		elapsed = unname(elapsed), output_value = store$value)
}

rd_dontrun_topics <- function() {
	rd_files <- sort(Sys.glob(file.path(root, "man", "*.Rd")))
	topics <- character()
	for(rd in rd_files) {
		lines <- readLines(rd, warn = FALSE)
		if(!any(grepl("\\\\dontrun\\{", lines))) next
		name_line <- grep("^\\\\name\\{", lines, value = TRUE)
		if(length(name_line) == 0) next
		topic <- sub("^\\\\name\\{([^}]+)\\}.*$", "\\1", name_line[1])
		topics <- c(topics, topic)
	}
	unique(topics)
}

run_example_with_dontrun <- function(topic) {
	env <- new.env(parent = globalenv())
	result <- capture_run(utils::example(topic, package = "depmixS4",
		ask = FALSE, echo = FALSE, local = env, run.dontrun = TRUE,
		character.only = TRUE))
	metrics <- collect_env_metrics(env)
	metrics$run_dontrun <- TRUE
	list(status = if(is.null(result$error)) "OK" else "FAIL",
		elapsed = result$elapsed, output = result$output,
		warnings = result$warnings, error = result$error, metrics = metrics)
}

dontrun_skip_reasons <- c(
	sp500 = "external TTR/Yahoo data-generation example"
)

dontrun_skip_reason <- function(topic) {
	if(topic %in% names(dontrun_skip_reasons)) {
		return(dontrun_skip_reasons[[topic]])
	}
	if(identical(topic, "makeDepmix") &&
		(!requireNamespace("gamlss", quietly = TRUE) ||
			!requireNamespace("gamlss.dist", quietly = TRUE))) {
		return("optional gamlss packages unavailable")
	}
	NULL
}

run_all <- function() {
	items <- empty_items()
	metrics <- empty_metrics()

	add_result <- function(suite, item, result) {
		out_name <- paste0(safe_name(paste(suite, item, sep = "-")), ".txt")
		out_path <- file.path(output_dir, out_name)
		error <- result$error
		if(is.null(error)) error <- ""
		output <- c(result$output,
			if(length(result$warnings)) c("", "Warnings:", result$warnings) else character(),
			if(nzchar(error)) c("", "Error:", error) else character())
		write_text(out_path, output)
		items <<- rbind(items, data.frame(suite = suite, item = item,
			status = result$status, elapsed_sec = result$elapsed,
			output_file = out_path, error = error,
			stringsAsFactors = FALSE))
		metrics <<- rbind(metrics,
			metric_rows(suite, item, result$metrics))
	}

	test_files <- sort(Sys.glob(file.path(root, "tests", "*.R")))
	for(test in test_files) {
		add_result("test", basename(test), run_test_file(test))
	}

	for(topic in rd_example_topics()) {
		add_result("example", topic, run_example(topic))
	}

	for(topic in rd_dontrun_topics()) {
		reason <- dontrun_skip_reason(topic)
		if(!is.null(reason)) {
			add_result("dontrun", topic,
				list(status = "SKIP", elapsed = 0, output = character(),
					warnings = character(), error = reason, metrics = list()))
		} else {
			add_result("dontrun", topic, run_example_with_dontrun(topic))
		}
	}

	list(items = items, metrics = metrics)
}

csv_path <- function(dir, name) file.path(dir, name)

write_results <- function(dir, data) {
	dir.create(dir, recursive = TRUE, showWarnings = FALSE)
	write.csv(data$items, csv_path(dir, "item_results.csv"), row.names = FALSE)
	write.csv(data$metrics, csv_path(dir, "metric_results.csv"), row.names = FALSE)
}

read_results <- function(dir) {
	list(
		items = read.csv(csv_path(dir, "item_results.csv"),
			stringsAsFactors = FALSE),
		metrics = read.csv(csv_path(dir, "metric_results.csv"),
			stringsAsFactors = FALSE)
	)
}

compare_metrics <- function(current, baseline, suite, item) {
	cur <- current[current$suite == suite & current$item == item, , drop = FALSE]
	base <- baseline[baseline$suite == suite & baseline$item == item, , drop = FALSE]
	keys <- sort(unique(c(cur$metric, base$metric)))
	if(length(keys) == 0) return("none")
	deviations <- character()
	for(key in keys) {
		crow <- cur[cur$metric == key, , drop = FALSE]
		brow <- base[base$metric == key, , drop = FALSE]
		if(nrow(crow) == 0) {
			deviations <- c(deviations, paste0(key, " missing current"))
		} else if(nrow(brow) == 0) {
			deviations <- c(deviations, paste0(key, " new"))
		} else if(!is.na(crow$numeric_value[1]) && !is.na(brow$numeric_value[1])) {
			diff <- crow$numeric_value[1] - brow$numeric_value[1]
			if(abs(diff) > 1e-6) {
				deviations <- c(deviations,
					paste0(key, " diff ", sprintf("%.6g", diff)))
			}
		} else if(!identical(crow$value[1], brow$value[1])) {
			deviations <- c(deviations, paste0(key, " changed"))
		}
	}
	if(length(deviations) == 0) "none" else paste(deviations, collapse = "; ")
}

make_comparison <- function(current, baseline) {
	cur <- current$items
	base <- baseline$items
	keys <- unique(rbind(cur[, c("suite", "item")], base[, c("suite", "item")]))
	rows <- lapply(seq_len(nrow(keys)), function(i) {
		suite <- keys$suite[i]
		item <- keys$item[i]
		crow <- cur[cur$suite == suite & cur$item == item, , drop = FALSE]
		brow <- base[base$suite == suite & base$item == item, , drop = FALSE]
		current_status <- if(nrow(crow)) crow$status[1] else "MISSING"
		baseline_status <- if(nrow(brow)) brow$status[1] else "MISSING"
		current_elapsed <- if(nrow(crow)) crow$elapsed_sec[1] else NA_real_
		baseline_elapsed <- if(nrow(brow)) brow$elapsed_sec[1] else NA_real_
		status_dev <- if(identical(current_status, baseline_status)) {
			"none"
		} else {
			paste0("status ", baseline_status, " -> ", current_status)
		}
		metric_dev <- compare_metrics(current$metrics, baseline$metrics,
			suite, item)
		deviation <- paste(c(status_dev, metric_dev)[c(status_dev, metric_dev) != "none"],
			collapse = "; ")
		if(!nzchar(deviation)) deviation <- "none"
		data.frame(
			suite = suite,
			item = item,
			current_status = current_status,
			baseline_status = baseline_status,
			current_elapsed_sec = current_elapsed,
			baseline_elapsed_sec = baseline_elapsed,
			elapsed_delta_sec = current_elapsed - baseline_elapsed,
			speed_ratio = current_elapsed / baseline_elapsed,
			runtime_improvement_pct = 100 * (baseline_elapsed - current_elapsed) /
				baseline_elapsed,
			deviation = deviation,
			stringsAsFactors = FALSE
		)
	})
	do.call(rbind, rows)
}

md_escape <- function(x) {
	x <- as.character(x)
	x <- gsub("\\|", "\\\\|", x)
	x <- gsub("\n", " ", x)
	x
}

format_num <- function(x) {
	ifelse(is.na(x), "", sprintf("%.3f", x))
}

write_markdown <- function(path, comparison, current_total, baseline_total,
	initialized_baseline) {
	lines <- c(
		"# depmixS4 Refactor Check Results",
		"",
		paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
		"",
		"Positive improvement percentages mean the current run is faster than the baseline.",
		"Refactor-check timings are single-run regression-check timings; use them as coarse signals, not benchmark-grade estimates.",
		"",
		"| Suite | Item | Current | Baseline | Current sec | Baseline sec | Delta sec | Speed ratio | Improvement % | Deviation |",
		"|---|---|---:|---:|---:|---:|---:|---:|---:|---|"
	)
	for(i in seq_len(nrow(comparison))) {
		row <- comparison[i, ]
		lines <- c(lines, paste(
			"|", md_escape(row$suite),
			"|", md_escape(row$item),
			"|", md_escape(row$current_status),
			"|", md_escape(row$baseline_status),
			"|", format_num(row$current_elapsed_sec),
			"|", format_num(row$baseline_elapsed_sec),
			"|", format_num(row$elapsed_delta_sec),
			"|", format_num(row$speed_ratio),
			"|", format_num(row$runtime_improvement_pct),
			"|", md_escape(row$deviation), "|"
		))
	}
	note <- paste0("Note: current total runtime was ",
		sprintf("%.3f", current_total), " seconds; baseline total runtime was ",
		sprintf("%.3f", baseline_total), " seconds.")
	note <- paste0(note,
		" Dontrun rows are generated from Rd topics containing ",
		"\\dontrun{} and run with utils::example(..., run.dontrun = TRUE).")
	if(initialized_baseline) {
		note <- paste0(note,
			" Baseline was initialized from this run; use --update-baseline ",
			"only when intentional behavior and speed changes should become ",
			"the new reference.")
	}
	lines <- c(lines, "", note)
	write_text(path, lines)
}

write_html <- function(md_path, html_path) {
	if(!requireNamespace("markdown", quietly = TRUE)) return(invisible(FALSE))
	body <- markdown::mark(file = md_path, format = "html", template = FALSE)
	write_text(html_path, c(
		"<!doctype html>",
		"<html>",
		"<head>",
		"<meta charset=\"utf-8\">",
		"<title>depmixS4 Refactor Check Results</title>",
		"<style>",
		"body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:2rem;line-height:1.45;}",
		"table{border-collapse:collapse;margin:1rem 0;}",
		"th,td{border:1px solid #ddd;padding:0.35rem 0.5rem;}",
		"th{background:#f6f8fa;}",
		"</style>",
		"</head>",
		"<body>",
		body,
		"</body>",
		"</html>"
	))
	invisible(TRUE)
}

message("Running refactor checks against depmixS4 ",
	utils::packageVersion("depmixS4"))
current <- run_all()
write_results(run_dir, current)
write_results(current_dir, current)

baseline_exists <- file.exists(csv_path(baseline_dir, "item_results.csv")) &&
	file.exists(csv_path(baseline_dir, "metric_results.csv"))
initialized_baseline <- FALSE
if(update_baseline || !baseline_exists) {
	write_results(baseline_dir, current)
	initialized_baseline <- !baseline_exists
}
baseline <- read_results(baseline_dir)
comparison <- make_comparison(current, baseline)
write.csv(comparison, file.path(current_dir, "comparison_table.csv"),
	row.names = FALSE)
write.csv(comparison, file.path(run_dir, "comparison_table.csv"),
	row.names = FALSE)

current_total <- sum(current$items$elapsed_sec, na.rm = TRUE)
baseline_total <- sum(baseline$items$elapsed_sec, na.rm = TRUE)
write_markdown(file.path(current_dir, "summary.md"), comparison,
	current_total, baseline_total, initialized_baseline)
write_markdown(file.path(run_dir, "summary.md"), comparison,
	current_total, baseline_total, initialized_baseline)
write_html(file.path(current_dir, "summary.md"),
	file.path(current_dir, "summary.html"))
write_html(file.path(run_dir, "summary.md"),
	file.path(run_dir, "summary.html"))

message("Wrote summary table to: ", file.path(current_dir, "summary.md"))
message("Current total runtime: ", sprintf("%.3f", current_total), " seconds")
message("Baseline total runtime: ", sprintf("%.3f", baseline_total), " seconds")
if(any(comparison$current_status != "OK" & comparison$current_status != "SKIP") ||
	any(comparison$deviation != "none")) {
	quit(status = 1)
}
