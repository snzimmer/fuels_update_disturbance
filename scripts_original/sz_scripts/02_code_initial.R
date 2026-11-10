library(terra)
library(data.table)
#
setwd(this.path::here())



suggest_rules <- function(status, rules) {
  
  # 1. Cell-level summary: leftover pixels, and whether a base rule already covers the cell
  cells <- status[, .(pixels = sum(N[status == "left_behind"]),
                      taken  = any(!is.na(top_spec) & top_spec == 0)),
                  by = .(fvt, fdist, fvc, fvh)]
  cells[, lf := paste(fvc %/% 10, fvh %/% 100)]
  cells <- cells[, if (sum(pixels) > 0) .SD, by = .(fvt, fdist, lf)]
  
  # 2. Greedy rectangle growth within one group
  grow <- function(g) {
    cv <- sort(unique(g$fvc)); hv <- sort(unique(g$fvh))
    M  <- matrix(0, length(cv), length(hv))
    Tk <- matrix(FALSE, length(cv), length(hv))
    idx <- cbind(match(g$fvc, cv), match(g$fvh, hv))
    M[idx] <- g$pixels; Tk[idx] <- g$taken
    out <- list()
    while (any(M > 0)) {
      s  <- which(M == max(M), arr.ind = TRUE)[1, ]
      r1 <- r2 <- s[1]; c1 <- c2 <- s[2]
      repeat {
        cand <- list(c(r1 - 1, r2, c1, c2), c(r1, r2 + 1, c1, c2),
                     c(r1, r2, c1 - 1, c2), c(r1, r2, c1, c2 + 1))
        gain <- sapply(cand, function(z) {
          if (z[1] < 1 || z[2] > nrow(M) || z[3] < 1 || z[4] > ncol(M)) return(-1)
          if (any(Tk[z[1]:z[2], z[3]:z[4]])) return(-1)
          sum(M[z[1]:z[2], z[3]:z[4]])
        })
        cur <- sum(M[r1:r2, c1:c2])
        if (max(gain) <= cur) break
        z <- cand[[which.max(gain)]]
        r1 <- z[1]; r2 <- z[2]; c1 <- z[3]; c2 <- z[4]
      }
      out[[length(out) + 1]] <- data.table(
        Cover_Low = cv[r1], Cover_High = cv[r2],
        Height_Low = hv[c1], Height_High = hv[c2],
        pixels = sum(M[r1:r2, c1:c2]))
      M[r1:r2, c1:c2] <- 0
    }
    rbindlist(out)
  }
  
  sug <- cells[, grow(.SD), by = .(fvt, fdist, lf)]
  
  # 3. Borrow fuel settings from the nearest existing rule (same FVT; prefer same FDist)
  hr  <- sort(unique(c(cells$fvh, rules$Height_Low, rules$Height_High)))
  hrk <- function(x) match(x, hr)
  act <- rules[on == TRUE]
  
  pick <- function(i) {
    s <- sug[i]
    d <- act[EVT == s$fvt & DIST %in% c(s$fdist, 0) &
               Cover_Low %/% 10 == s$Cover_Low %/% 10]
    if (!nrow(d)) d <- act[EVT == s$fvt]
    if (!nrow(d)) return(data.table(donor_ID = NA_integer_))
    dist <- abs((d$Cover_Low + d$Cover_High) / 2 - (s$Cover_Low + s$Cover_High) / 2) +
      abs((hrk(d$Height_Low) + hrk(d$Height_High)) / 2 -
            (hrk(s$Height_Low) + hrk(s$Height_High)) / 2) +
      100 * (d$DIST != s$fdist)
    k <- d[which.min(dist)]
    data.table(donor_ID = k$ID, FBFM13 = k$FBFM13, FBFM40 = k$FBFM40,
               Canopy = k$Canopy, FCCS = k$FCCS)
  }
  donors <- rbindlist(lapply(seq_len(nrow(sug)), pick), fill = TRUE)
  
  draft <- cbind(sug, donors)
  draft[, `:=`(EVT = fvt, DIST = fdist, BPSRF = "any", Wildcard = "any",
               OnOff = "On",
               Notes = paste0("Draft rule; settings copied from rule ", donor_ID))]
  setorder(draft, fvt)
  draft[, .(EVT, DIST, Cover_Low, Cover_High, Height_Low, Height_High,
            BPSRF, Wildcard, FBFM13, FBFM40, Canopy, OnOff, Notes,
            pixels, donor_ID)]
}

# Combine input rasters and get pixel counts
fvt<- rast("./inputs/Chloe/OkaWen_LFTFC/LFTFC_MaxEffect/Input/Zone_01/FOA/landfire2020_FOA_Z01_FVT_120m.tif")
fvc<- rast("./inputs/Chloe/OkaWen_LFTFC/LFTFC_MaxEffect/Input/Zone_01/FOA/landfire2020_FOA_Z01_FVC_120m.tif")
fvh<- rast("./inputs/Chloe/OkaWen_LFTFC/LFTFC_MaxEffect/Input/Zone_01/FOA/landfire2020_FOA_Z01_FVH_120m.tif")
fdist<- rast("./inputs/Chloe/OkaWen_LFTFC/LFTFC_MaxEffect/Input/Zone_01/FOA/landfire2020_FOA_Z01_FDist121_120m.tif")
bps<- rast("./inputs/Chloe/OkaWen_LFTFC/LFTFC_MaxEffect/Input/Zone_01/FOA/landfire2020_FOA_Z01_BPS_120m.tif")


stack<- c(fvt, fvc, fvh, fdist, bps)      
names(stack)<- c("fvt","fvc","fvh","fdist", "bps")

# strip levels
for (i in 1:nlyr(stack)) levels(stack[[i]]) <- NULL 


count_combos <- function(r) {
  b <- terra::blocks(r)
  readStart(r); on.exit(readStop(r))
  
  parts <- vector("list", b$n)
  for (i in seq_len(b$n)) {
    v <- as.data.table(readValues(r, row = b$row[i], nrows = b$nrows[i],
                                  dataframe = TRUE))
    v <- v[complete.cases(v)]
    parts[[i]] <- v[, .N, by = names(v)]
  }
  rbindlist(parts)[, .(N = sum(N)), by = names(r)]
}

combos <- count_combos(stack)
combos[, combo_id := .I]



# Read in rule for the zone
rules<- as.data.table(read.csv("./LFTFC/tables/z01_nw24_250_Rulesets.csv"))
# remove undisturbed rules
rules<- rules[rules$DIST != 0 ,]

#
# Rename columns
rules[, `:=`(
  fvt        = EVT,
  dist       = DIST,
  bps        = fifelse(BPSRF    == "any", NA_character_, as.character(BPSRF)),
  wild       = fifelse(Wildcard == "any", NA_character_, as.character(Wildcard)),
  on         = OnOff == "On",
  fbfm40_lab = trimws(sub("/.*", "", FBFM40)),
  fbfm40     = as.integer(sub(".*/\\s*", "", FBFM40))
)]


# Columns that are plain numbers with a 9999 sentinel
num_nodata <- c("CCover","CHeight","CBD13x100","CBD40x100","CBH13mx10","CBH40mx10", "FLM")
rules[, (num_nodata) := lapply(.SD, function(x) {
  x <- as.numeric(x)
  fifelse(x == 9999, NA_real_, x)
}), .SDcols = num_nodata]

# Label columns in "LABEL / code" form (FBFM40, CanFM)
split_label <- function(x) {
  x <- as.character(x)
  x[x == "9999"] <- NA
  list(label = trimws(sub("\\s*/.*", "", x)),
       code  = suppressWarnings(as.integer(sub(".*/\\s*", "", x))))
}

fb <- split_label(rules$FBFM40)
rules[, `:=`(fbfm40_lab = fb$label, fbfm40 = fb$code)]

cf <- split_label(rules$CanFM)
rules[, `:=`(canfm_lab = cf$label, canfm = cf$code)]
#
#rules$bps<- NA
#rules$wild<- NA

# Specificity per guide 6.5.3: base < BPS-only / wildcard-only < both
rules[, spec := fifelse(!is.na(bps) & !is.na(wild), 2L,
                        fifelse(!is.na(bps) | !is.na(wild), 1L, 0L))]



#rules<- rules[-which(rules$ID %in% c(8656:8700, 1:500)),]


rr <- rules[on == TRUE,
            .(rule_id = ID, fvt, dist, Cover_Low, Cover_High,
              Height_Low, Height_High, spec,
              rbps  = suppressWarnings(as.integer(bps)),
              rwild = suppressWarnings(as.integer(wild)))]
if (!"wild" %in% names(combos)) combos[, wild := NA_integer_]


# 3. Match: each combo to every rule whose FVT, DIST, cover and height fit
m <- rr[combos,
        on = .(fvt, dist = fdist,
               Cover_Low  <= fvc, Cover_High  >= fvc,
               Height_Low <= fvh, Height_High >= fvh),
        nomatch = NULL, allow.cartesian = TRUE]

# 4. BPS / wildcard restrictions (bps and wild here are the combos' values)
m <- m[(is.na(rbps) | rbps == bps) & (is.na(rwild) | rwild == wild)]

# 5. Pick the most specific rule per combo; flag ties as overlaps
best <- m[, {
  top <- max(spec)
  .(top_spec = top,
    n_top    = sum(spec == top),
    rule_id  = rule_id[which.max(spec)])
}, by = combo_id]

status <- merge(combos, best, by = "combo_id", all.x = TRUE)
status[, status := fcase(is.na(top_spec), "left_behind",
                         n_top > 1,       "overlap",
                         default =        "ok")]

# Summary (the tool's red/green label, but for the whole MU)
status[, .(combos = .N, pixels = sum(N)), by = status]

# Assign FBFM for the combos with existing rules
status$fbfm<-NA
status$fbfm_lab<-NA
#
status$fbfm[status$status == "ok"]<- rules$fbfm40[match(status$rule_id[status$status == "ok"], rules$ID)]
status$fbfm_lab[status$status == "ok"]<- rules$fbfm40_lab[match(status$rule_id[status$status == "ok"], rules$ID)]



# Find left-behind pixels needing rules
left_behind<- status[status == "left_behind",
                     .(pixels = sum(N)),
                     by = .(fvt, fvc, fvh, fdist, bps)][order(fvt)]




potential_rules <- suggest_rules(status, rules)
potential_rules
