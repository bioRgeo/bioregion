exportGDF <- function(df,
                      col1 = "Node1",
                      col2 = "Node2",
                      weight = NULL,
                      node_chars = NULL,
                      node_chars_ID = colnames(node_chars)[1],
                      color_column = NULL,
                      file = "output.gdf") {
  
  nodes <- unique(c(df[[col1]], df[[col2]]))

  node_df <- data.frame(name = nodes, label = nodes, stringsAsFactors = FALSE)

  if (!is.null(node_chars)) {
    colnames(node_chars)[which(colnames(node_chars) ==
                                 node_chars_ID)] <- "name"
    if (!('name' %in% colnames(node_chars))) {
      stop("The node_chars data.frame must contain a 'name' column.")
    }
     node_df <- merge(node_df, node_chars, by = 'name', all.x = TRUE)
  }
  
  if (!is.null(color_column)) {
    if (!color_column %in% names(node_df)) {
      stop(paste("Color column", color_column, "not found in node characteristics"))
    }
    
    # Convert the color column from hex to RGB
    hex_to_rgb <- function(hex) {
      # Remove '#' if present
      hex <- gsub("#", "", hex)
      # Split into R, G, B components
      rgb_values <- sapply(seq(1, nchar(hex), 2), function(i) substr(hex, i, i+1))
      # Convert hex to decimal
      rgb_values <- as.numeric(sapply(rgb_values, function(x) strtoi(x, 16L)))
      # Return as a string "R,G,B"
      paste(rgb_values, collapse=",")
    }
    
    node_df$color <- sapply(node_df[[color_column]], hex_to_rgb)
    # Remove the original color column if not needed
    node_df[[color_column]] <- NULL
  }
  
  node_attributes <- setdiff(names(node_df), c("name", "label"))
  attribute_definitions <- sapply(node_attributes, function(attr) {
    if (attr == "color") {
      "graphics.color VARCHAR"
    } else if (is.numeric(node_df[[attr]])) {
      paste0(attr, " DOUBLE")
    } else {
      paste0(attr, " VARCHAR")
    }
  })
  
  node_header <- paste("nodedef>name VARCHAR,label VARCHAR", paste(attribute_definitions, collapse=","), sep=",")
  
  node_lines <- apply(node_df, 1, function(x) paste(x, collapse=","))
  
  if (!is.null(weight)) {
    edge_header <- "edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE"
    edge_lines <- apply(df[, c(col1, col2, weight)], 1, paste, collapse=",")
  } else {
    edge_header <- "edgedef>node1 VARCHAR,node2 VARCHAR"
    edge_lines <- apply(df[, c(col1, col2)], 1, paste, collapse=",")
  }
  
  all_lines <- c(node_header, node_lines, edge_header, edge_lines)
  writeLines(all_lines, con = file)
}
