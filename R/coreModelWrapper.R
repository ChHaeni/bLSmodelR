coreModelWrapper <- function(uvw, sncrun) {
	coreModel(uvw[, 'u0'], uvw[, 'v0'], uvw[, 'w0'], sncrun[, SensorHeight], 
        sncrun[, Ustar], sncrun[, L], sncrun[, Zo], sncrun[, bw], sncrun[, sUu], sncrun[, sVu], sncrun[, kv], sncrun[, C0], sncrun[, alpha], sncrun[, MaxFetch])
}
