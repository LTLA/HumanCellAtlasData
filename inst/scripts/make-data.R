library(rhdf5)
library(R.utils)

tsv2hdf5 <- function(path, output.file, output.name) {
    # Getting the number of lines (minus header)
    nlines <- countLines(path)
    nlines <- as.integer(nlines - 1L)

    # Getting the header and the number of fields.
    incon <- file(path)
    open(incon)
    on.exit({ close(incon) })

    header <- readLines(incon, 1)
    col.names <- strsplit(header, "\t")[[1]][-1] # getting rid of empty row name.
    nfields <- length(col.names)

    # Constructing the output HDF5 file.
    h5createFile(output.file)
    chunking <- beachmat::getBestChunkDims(c(nlines, nfields))
    h5createDataset(output.file, output.name, 
                    dims=c(nlines, nfields), 
                    storage.mode="integer",
                    chunk=chunking)

    h5con <- H5Fopen(output.file)
    on.exit({ H5Fclose(h5con) }, add=TRUE)

    # Looping over the input TSV file and pulling out chunks for writing.
    cur.row <- 1L
    row.names <- vector("list", ceiling(nlines/chunking[1]))
    col.what <- c(list("character"), as.list(integer(nfields)))
    idx <- 1L

    while (cur.row <= nlines) { 
        current.lot <- scan(incon, nlines=chunking[1], what=col.what) 
        row.names[[idx]] <- current.lot[[1]]
        idx <- idx + 1L

        current.lot <- current.lot[-1]
        current.lot <- do.call(cbind, current.lot)
        h5writeDataset(current.lot, h5con, output.name, start=c(cur.row, 1L))
        cur.row <- cur.row + nrow(current.lot)
    }
  
    # Returning dimnames.
    return(list(rownames=unlist(row.names), colnames=col.names))
}

# # Creating a dummy file.
# blah <- matrix(1:100000, nrow=100)
# rownames(blah) <- paste0("Gene", 1:nrow(blah))
# colnames(blah) <- paste0("Cell", 1:ncol(blah))
# rownames(blah) <- paste0("Gene", 1:nrow(blah))
# colnames(blah) <- paste0("Cell", 1:ncol(blah))
# output <- gzfile("whee.tsv.gz", open="wb")
# write.table(blah, file=output, sep="\t", col.names=NA, quote=FALSE)
# close(output)
#
# # Using the above function to store it in a HDF5 file.
# system.time(blah <- tsv2hdf5("whee.tsv.gz", "whee.h5", "counts"))
# library(HDF5Array)
# h5mat <- HDF5Array("whee.h5", "counts")
# dimnames(h5mat) <- blah

library(HDF5Array)
library(SingleCellExperiment)
loom2sce <- function(path) {
    # Extracting the count matrix.
    mat <- HDF5Array(path, "matrix")
    mat <- t(mat)

    # Extracting the row and column metadata.
    col.attrs <- h5read(path, "col_attrs")
    if (length(col.attrs)) { 
        col.df <- data.frame(col.attrs)
    } else {
        col.df <- DataFrame(matrix(0, nrow(mat), 0))
    }

    row.attrs <- h5read(path, "row_attrs")
    if (length(row.attrs)) { 
        row.df <- data.frame(row.attrs)
    } else {
        row.df <- NULL
    }

    # Extracting layers (if there are any)
    optional <- h5ls(path)
    is.layer <- optional$group=="/layer"
    if (any(is.layer)) {
        layer.names <- optional$name[is.layer,]
        other.layers <- vector("list", length(layer.names))
        names(other.layers) <- layer.names

        for (layer in layer.names) {
            current <- HDF5Array(path, file.path("/layer", layer))
            other.layers[[layer]] <- t(current)
        }
    } else {
        other.layers <- list()
    }
    
    # Returning SingleCellExperiment object.
    sce <- SingleCellExperiment(c(matrix=mat, other.layers), rowData=row.df, colData=col.df)
    return(sce)
}

# sce <- loom2sce("cortex.loom")
