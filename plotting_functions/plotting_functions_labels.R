prepare_label_contours <- function(file, slice_axis, slice_axis_coordinates, labels=NULL, mask_file=NULL, filter_size=NULL) {

  # Get file shape
  file_shape <- get_file_shape(file)

  # Get voxel indices
  slices_df <- get_slice_mappings(file = file, slice_axis = slice_axis, slice_axis_coordinates = slice_axis_coordinates)

  # Read in 3D file
  arr <- round(mincArray(mincGetVolume(file)))
  if (!is.null(mask_file)) {
    maskarr <- mincArray(mincGetVolume(mask_file))
  }

  # Get array squence based on step direction
  # This is so that we can flip array based on direction
  # (a workaround since contourLines requires x, y array coordinates be provided in ascending order, but step direction may be negative)
  axis_sequence <- 1:3 %>%
    map(function(ax) {
      ax_dir <- file_shape$direction[[ax]]
      if (ax_dir == -1) {
        ax_seq <- rev(1:dim(arr)[[ax]])
      } else {
        ax_seq <- 1:dim(arr)[[ax]]
      }
    })

  # Get slices
  contours_df <- seq_along(slices_df$slice_index) %>%
    map_dfr(
      function(i) {

        # Get world and voxel coordinates corresponding to slice
        w <- slices_df$world[i]
        v <- slices_df$voxel[i]
        index <- slices_df$slice_index[i]

        # Get 2D array corresponding to background intensities
        switch (slice_axis,
                x={

                  # Get slice array
                  slice_array <- arr[(v + 1),,]
                  if (!is.null(mask_file)) {
                    mask_array <- maskarr[(v + 1),,]
                    slice_array[mask_array < 0.5] <- 0
                  }

                  # Flip array according to direction
                  slice_array <- slice_array[axis_sequence[[2]], axis_sequence[[3]]]
                  if (!is.null(mask_file)) {
                    mask_array <- mask_array[axis_sequence[[2]], axis_sequence[[3]]]
                  }

                  # Get unique label values
                  slice_label_values <- sort(unique(matrix(slice_array)))
                  if (!is.null(labels)) {
                    slice_label_values <- intersect(slice_label_values, labels)
                  }
                  slice_label_values <- setdiff(slice_label_values, 0)

                  # Loop over labels and get contours
                  out <- slice_label_values %>%
                    map_dfr(function(lv) {
                      slice_label_array <- slice_array
                      slice_label_array[] <- 0
                      slice_label_array[slice_array==lv] <- 1

                      # Use grDevices::contourLines() to get the contours
                      # Returns a list of contour paths
                      slice_contours <- contourLines(x = seq(from=file_shape$starts["y"],
                                                             to=file_shape$ends["y"],
                                                             length.out = dim(arr)[[2]]),
                                                     y = seq(from=file_shape$starts["z"],
                                                             to=file_shape$ends["z"],
                                                             length.out = dim(arr)[[3]]),
                                                     z = slice_label_array,
                                                     levels = 0.5)

                      # Wrangle data to long form
                      outlabel <- 1:length(slice_contours) %>%
                        map(function(i) {
                          as_tibble(slice_contours[[i]]) %>%
                            mutate(obj=i)
                        }
                        ) %>%
                        bind_rows() %>%
                        rename(y=x, z=y) %>%
                        mutate(x=w,
                               slice_axis=slice_axis,
                               slice_voxel=v,
                               slice_world=w,
                               slice_index=index,
                               label=lv) %>%
                        select(label, obj, x, y, z, slice_axis, slice_voxel, slice_world, slice_index)

                      return(outlabel)
                    })
                },
                y={

                  # Get slice array
                  slice_array <- arr[,(v + 1),]
                  if (!is.null(mask_file)) {
                    mask_array <- maskarr[,(v + 1),]
                    slice_array[mask_array < 0.5] <- 0
                  }

                  # Flip array according to direction
                  slice_array <- slice_array[axis_sequence[[1]], axis_sequence[[3]]]
                  if (!is.null(mask_file)) {
                    mask_array <- mask_array[axis_sequence[[1]], axis_sequence[[3]]]
                  }

                  # Get unique label values
                  slice_label_values <- sort(unique(matrix(slice_array)))
                  if (!is.null(labels)) {
                    slice_label_values <- intersect(slice_label_values, labels)
                  }
                  slice_label_values <- setdiff(slice_label_values, 0)

                  # Loop over labels and get contours
                  out <- slice_label_values %>%
                    map_dfr(function(lv) {
                      slice_label_array <- slice_array
                      slice_label_array[] <- 0
                      slice_label_array[slice_array==lv] <- 1

                      # Use grDevices::contourLines() to get the contours
                      # Returns a list of contour paths
                      slice_contours <- contourLines(x = seq(from=file_shape$starts["x"],
                                                             to=file_shape$ends["x"],
                                                             length.out = dim(arr)[[1]]),
                                                     y = seq(from=file_shape$starts["z"],
                                                             to=file_shape$ends["z"],
                                                             length.out = dim(arr)[[3]]),
                                                     z = slice_label_array,
                                                     levels = 0.5)

                      # Wrangle data to long form
                      outlabel <- 1:length(slice_contours) %>%
                        map(function(i) {
                          as_tibble(slice_contours[[i]]) %>%
                            mutate(obj=i)
                        }
                        ) %>%
                        bind_rows() %>%
                        rename(x=x, z=y) %>%
                        mutate(y=w,
                               slice_axis=slice_axis,
                               slice_voxel=v,
                               slice_world=w,
                               slice_index=index,
                               label=lv) %>%
                        select(label, obj, x, y, z, slice_axis, slice_voxel, slice_world, slice_index)

                      return(outlabel)
                    })



                },
                z={

                  # Get slice array
                  slice_array <- arr[,,(v + 1)]
                  if (!is.null(mask_file)) {
                    mask_array <- maskarr[,,(v + 1)]
                    slice_array[mask_array < 0.5] <- 0
                  }

                  # Flip array according to direction
                  slice_array <- slice_array[axis_sequence[[1]], axis_sequence[[2]]]
                  if (!is.null(mask_file)) {
                    mask_array <- mask_array[axis_sequence[[1]], axis_sequence[[2]]]
                  }

                  # Get unique label values
                  slice_label_values <- sort(unique(matrix(slice_array)))
                  if (!is.null(labels)) {
                    slice_label_values <- intersect(slice_label_values, labels)
                  }
                  slice_label_values <- setdiff(slice_label_values, 0)

                  # Loop over labels and get contours
                  out <- slice_label_values %>%
                    map_dfr(function(lv) {
                      slice_label_array <- slice_array
                      slice_label_array[] <- 0
                      slice_label_array[slice_array==lv] <- 1

                      # Use grDevices::contourLines() to get the contours
                      # Returns a list of contour paths
                      slice_contours <- contourLines(x = seq(from=file_shape$starts["x"],
                                                             to=file_shape$ends["x"],
                                                             length.out = dim(arr)[[1]]),
                                                     y = seq(from=file_shape$starts["y"],
                                                             to=file_shape$ends["y"],
                                                             length.out = dim(arr)[[2]]),
                                                     z = slice_label_array,
                                                     levels = 0.5)

                      # Wrangle data to long form
                      outlabel <- 1:length(slice_contours) %>%
                        map(function(i) {
                          as_tibble(slice_contours[[i]]) %>%
                            mutate(obj=i)
                        }
                        ) %>%
                        bind_rows() %>%
                        rename(x=x, y=y) %>%
                        mutate(z=w,
                               slice_axis=slice_axis,
                               slice_voxel=v,
                               slice_world=w,
                               slice_index=index,
                               label=lv) %>%
                        select(label, obj, x, y, z, slice_axis, slice_voxel, slice_world, slice_index)

                      return(outlabel)
                    })
                }
        )



        # Return
        return(out)
      }
    )

  # Filter by size if required
  # Many contour paths exist. Remove the small ones by filtering out paths with few points
  if (!is.null(filter_size)) {
    contours_df <- contours_df %>%
      group_by(obj, slice_index) %>%
      filter(length(label) >= filter_size) %>%
      ungroup()
  }

  # Return
  out <- list(file=file,
              contours_df=contours_df,
              slice_axis=slice_axis,
              slice_axis_coordinates=slice_axis_coordinates,
              labels=sort(unique(contours_df$label)),
              mask_file=mask_file,
              filter_size=filter_size)
}
