library(tidyverse)
library(glue)
library(RMINC)
library(scales)
library(patchwork)
library(ggnewscale)

#############################

get_file_shape <- function(file) {

  # Get dimensions of array
  file_dims <- minc.dimensions.sizes(file) # Returns lengths as z, y, x
  names(file_dims) <- c("z", "y", "x")

  # Get step size (resolution) of array
  file_steps <- minc.separation.sizes(file) # Returns steps as z, y, x
  file_step_signs <- sign(file_steps) # Also z, y, x
  names(file_steps) <- c("z", "y", "x")
  names(file_step_signs) <- c("z", "y", "x")

  # Get coordinate ranges of array
  file_starts <- mincConvertVoxelToWorld(file, 0, 0, 0) # Input z, y, x; returns as x, y, z
  file_ends <- mincConvertVoxelToWorld(file, file_dims[1]-1, file_dims[2]-1, file_dims[3]-1) # Input z, y, x; returns as x, y, z
  file_range <- file_ends - file_starts
  names(file_starts) <- c("x", "y", "z")
  names(file_ends) <- c("x", "y", "z")
  names(file_range) <- c("x", "y", "z")

  # Object to return
  out <- list(file=file,
              dims=rev(file_dims),
              steps=rev(file_steps)*rev(file_step_signs),
              direction=rev(file_step_signs),
              starts=c(ifelse(file_step_signs["x"] > 0, file_starts["x"], file_ends["x"]),
                       ifelse(file_step_signs["y"] > 0, file_starts["y"], file_ends["y"]),
                       ifelse(file_step_signs["z"] > 0, file_starts["z"], file_ends["z"])
              ),
              ends=c(ifelse(file_step_signs["x"] > 0, file_ends["x"], file_starts["x"]),
                     ifelse(file_step_signs["y"] > 0, file_ends["y"], file_starts["y"]),
                     ifelse(file_step_signs["z"] > 0, file_ends["z"], file_starts["z"])
              ),
              range=rev(file_step_signs)*file_range
  )

  # Return
  return(out)
}

#############################

get_grid_sequence <- function(file, grid_padding=-0.5, line_points=100, grid_spacing=NULL, grid_lines=NULL) {

  file_shape <- get_file_shape(file)

  # Get starting coordinate
  start_x <- file_shape$starts["x"] + grid_padding
  start_y <- file_shape$starts["y"] + grid_padding
  start_z <- file_shape$starts["z"] + grid_padding


  # Get ending coordinate
  end_x <- file_shape$ends["x"] - grid_padding
  end_y <- file_shape$ends["y"] - grid_padding
  end_z <- file_shape$ends["z"] - grid_padding

  # Define position of bounding box for gridlines
  box_sequence <- list(x=c(start_x, end_x),
                       y=c(start_y, end_y),
                       z=c(start_z, end_z)
  )

  # Define position of gridlines along each axis
  if (!is.null(grid_lines) & !is.null(grid_spacing)) {
    stop("You must only specify exactly one of grid_lines or grid_spacing!")
  }

  if (!is.null(grid_lines) & is.null(grid_spacing)) {
    grid_sequence <- list(x=seq(from=start_x,
                                to=end_x,
                                length.out=grid_lines),
                          y=seq(from=start_y,
                                to=end_y,
                                length.out=grid_lines),
                          z=seq(from=start_z,
                                to=end_z,
                                length.out=grid_lines)
    )
  } else if (!is.null(grid_spacing) & is.null(grid_lines)) {
    grid_sequence <- list(x=seq(from=start_x,
                                to=end_x,
                                by=grid_spacing),
                          y=seq(from=start_y,
                                to=end_y,
                                by=grid_spacing),
                          z=seq(from=start_z,
                                to=end_z,
                                by=grid_spacing)
    )
  } else {
    stop("You must only specify exactly one of grid_lines or grid_spacing!")
  }


  # Define position of points along each grid and box line
  point_sequence <- list(x=seq(from=start_x,
                               to=end_x,
                               length.out=line_points),
                         y=seq(from=start_y,
                               to=end_y,
                               length.out=line_points),
                         z=seq(from=start_z,
                               to=end_z,
                               length.out=line_points)
  )

  # Output object
  out <- list(file=file,
              box_sequence=box_sequence,
              grid_sequence=grid_sequence,
              point_sequence=point_sequence)

  # Return
  return(out)
}

#############################

get_grid_sequence_manual <- function(file, starts, ends, line_points=100, grid_spacing=NULL, grid_lines=NULL) {

  file_shape <- get_file_shape(file)

  # Get starting coordinate
  start_x <- starts[[1]]
  start_y <- starts[[2]]
  start_z <- starts[[3]]


  # Get ending coordinate
  end_x <- ends[[1]]
  end_y <- ends[[2]]
  end_z <- ends[[3]]

  # Define position of bounding box for gridlines
  box_sequence <- list(x=c(start_x, end_x),
                       y=c(start_y, end_y),
                       z=c(start_z, end_z)
  )

  # Define position of gridlines along each axis
  if (!is.null(grid_lines) & !is.null(grid_spacing)) {
    stop("You must only specify exactly one of grid_lines or grid_spacing!")
  }

  if (!is.null(grid_lines) & is.null(grid_spacing)) {
    grid_sequence <- list(x=seq(from=start_x,
                                to=end_x,
                                length.out=grid_lines),
                          y=seq(from=start_y,
                                to=end_y,
                                length.out=grid_lines),
                          z=seq(from=start_z,
                                to=end_z,
                                length.out=grid_lines)
    )
  } else if (!is.null(grid_spacing) & is.null(grid_lines)) {
    grid_sequence <- list(x=seq(from=start_x,
                                to=end_x,
                                by=grid_spacing),
                          y=seq(from=start_y,
                                to=end_y,
                                by=grid_spacing),
                          z=seq(from=start_z,
                                to=end_z,
                                by=grid_spacing)
    )
  } else {
    stop("You must only specify exactly one of grid_lines or grid_spacing!")
  }

  # Define position of points along each grid and box line
  point_sequence <- list(x=seq(from=start_x,
                               to=end_x,
                               length.out=line_points),
                         y=seq(from=start_y,
                               to=end_y,
                               length.out=line_points),
                         z=seq(from=start_z,
                               to=end_z,
                               length.out=line_points)
  )

  # Output object
  out <- list(file=file,
              box_sequence=box_sequence,
              grid_sequence=grid_sequence,
              point_sequence=point_sequence)

  # Return
  return(out)
}

#############################

get_grid_sequence_around_point <- function(file, point, grid_spacing, grid_lines, line_points=100) {

  file_shape <- get_file_shape(file)

  # Get starts and ends
  grid_range <- grid_spacing*(grid_lines-1)

  start_x <- point[[1]] - grid_range/2
  start_y <- point[[2]] - grid_range/2
  start_z <- point[[3]] - grid_range/2

  end_x <- point[[1]] + grid_range/2
  end_y <- point[[2]] + grid_range/2
  end_z <- point[[3]] + grid_range/2

  out <- get_grid_sequence_manual(file=file,
                                  starts=c(start_x, start_y, start_z),
                                  ends=c(end_x, end_y, end_z),
                                  line_points=line_points, grid_lines=grid_lines)

  return(out)

}

#############################

# Function to create a series points along lines at a sequence of points
# Creates lines along one axis (half of what's needed for a rectangular grid array)
get_half_grid <- function(grid_coordinates, # Sequence of points along main "grid" axis
                          line_coordinates, # Sequence of points along each grid line
                          grid_type="grid", # Type of grid (full grid or boundingbox ?)
                          grid_space="native", # Space in which grid is defined
                          grid_axis=NA,  # Which axis is the "grid" axis?
                          slice_axis=NA, # Which axis is the slice axis?
                          slice_axis_coordinate=NA) {

  # Get axes
  all_axes <- c("x", "y", "z")
  line_axis <- setdiff(all_axes, c(grid_axis, slice_axis))

  # Map along each of the axis coordinates to draw a line
  out <- seq_along(grid_coordinates) %>%
    map_dfr(function(lc) {

      # Map over points in the line
      seq_along(line_coordinates) %>%
        map_dfr(function(i) {
          tibble(grid_type=grid_type,
                 grid_space=grid_space,
                 slice_axis=slice_axis,
                 grid_axis=grid_axis,
                 gridline=lc,
                 slice_axis_coordinate=slice_axis_coordinate,
                 grid_axis_coordinate=grid_coordinates[[lc]],
                 line_axis_coordinate=line_coordinates[[i]])
        })

    })

  # Rename axes
  out <- out %>%
    rename(!!sym(slice_axis):=slice_axis_coordinate,
           !!sym(grid_axis):=grid_axis_coordinate,
           !!sym(line_axis):=line_axis_coordinate) %>%
    select(grid_type, grid_space,
           slice_axis, grid_axis,
           gridline,
           x, y, z)

  # Return
  return(out)
}

#############################

get_full_grid <- function(grid_sequence,  # Sequence of points along main "grid" axis (for all axes)
                          point_sequence, # Sequence of points along each grid line (for all axes)
                          grid_type="grid", #  Type of grid (full grid or boundingbox ?) - depends on grid_sequence input type
                          grid_space=NA, # Space in which grid is defined)
                          slice_axis="y", # Slice plane
                          slice_axis_coordinate=NA # Which slice is this?
) {

  # Get axes
  all_axes <- c("x", "y", "z")
  grid_and_line_axes <- setdiff(all_axes, slice_axis)

  # Compute vertical grid
  grid_axis <- grid_and_line_axes[[1]]
  line_axis <- grid_and_line_axes[[2]]

  vertical_grid <- get_half_grid(grid_coordinates = grid_sequence[[grid_axis]],
                                 line_coordinates = point_sequence[[line_axis]],
                                 grid_type=grid_type,
                                 grid_space=grid_space,
                                 grid_axis=grid_axis,
                                 slice_axis=slice_axis,
                                 slice_axis_coordinate=slice_axis_coordinate)

  # Compute horizontal grid
  grid_axis <- grid_and_line_axes[[2]]
  line_axis <- grid_and_line_axes[[1]]

  horizontal_grid <- get_half_grid(grid_coordinates = grid_sequence[[grid_axis]],
                                   line_coordinates = point_sequence[[line_axis]],
                                   grid_type=grid_type,
                                   grid_space=grid_space,
                                   grid_axis=grid_axis,
                                   slice_axis=slice_axis,
                                   slice_axis_coordinate=slice_axis_coordinate)

  # Bind perpendicular grids together
  out <- rbind(vertical_grid, horizontal_grid)

  # Return
  return(out)

}

#############################

get_base_grid <- function(gs, # Output of get_grid_sequence()
                          grid_space=NA,
                          slice_axis="y",
                          slice_axis_coordinate=0
) {

  # Get the point sequences
  box_sequence <- gs$box_sequence
  grid_sequence <- gs$grid_sequence
  point_sequence <- gs$point_sequence

  # Get the grid and box sequence
  grid_df <- get_full_grid(grid_sequence = gs$grid_sequence,
                           point_sequence = gs$point_sequence,
                           grid_type = "grid",
                           grid_space = grid_space,
                           slice_axis = slice_axis,
                           slice_axis_coordinate = slice_axis_coordinate)

  box_df <- get_full_grid(grid_sequence = gs$box_sequence,
                          point_sequence = gs$point_sequence,
                          grid_type = "box",
                          grid_space = grid_space,
                          slice_axis = slice_axis,
                          slice_axis_coordinate = slice_axis_coordinate)

  # Output
  out <- list(grid=grid_df,
              box=box_df,
              file=gs$file,
              grid_space=grid_space,
              slice_axis=slice_axis,
              slice_axis_coordinate=slice_axis_coordinate,
              transformation_info=NA)

  # Return
  return(out)
}

#############################

get_slice_mappings <- function(file, slice_axis, slice_axis_coordinates) {
  slices_df <- slice_axis_coordinates %>%
    map_dfr(
      function(w) {
        world_coord <- w
        switch (slice_axis,
                x={
                  voxel_coord <- mincConvertWorldToVoxel(file, w, 0, 0)[3] # Input x,y,z -> returns z, y, x (zero indexed)
                },
                y={
                  voxel_coord <- mincConvertWorldToVoxel(file, 0, w, 0)[2] # Input x,y,z -> returns z, y, x (zero indexed)
                },
                z={
                  voxel_coord <- mincConvertWorldToVoxel(file, 0, 0, w)[1] # Input x,y,z -> returns z, y, x (zero indexed)
                }
        )

        out <- tibble(slice_axis=slice_axis,
                      world=world_coord,
                      voxel=voxel_coord)
        return(out)
      }
    ) %>%
    arrange(desc(world)) %>%
    mutate(slice_index=1:nrow(.)) %>%
    select(slice_axis, slice_index, world, voxel)

  # Return
  return(slices_df)
}


#############################

prepare_anatomy <- function(file, slice_axis, slice_axis_coordinates) {

  # Get file shape
  file_shape <- get_file_shape(file)

  # Get voxel indices
  slices_df <- get_slice_mappings(file = file, slice_axis = slice_axis, slice_axis_coordinates = slice_axis_coordinates)

  # Read in 3D file
  arr <- mincArray(mincGetVolume(file))

  # Set dimension names
  dimnames(arr) <- list(
    as.character(seq(from=ifelse(file_shape$direction["x"] > 0, file_shape$starts["x"], file_shape$ends["x"]),
                     to=ifelse(file_shape$direction["x"] > 0, file_shape$ends["x"], file_shape$starts["x"]),
                     length.out = dim(arr)[1])
    ),
    as.character(seq(from=ifelse(file_shape$direction["y"] > 0, file_shape$starts["y"], file_shape$ends["y"]),
                     to=ifelse(file_shape$direction["y"] > 0, file_shape$ends["y"], file_shape$starts["y"]),
                     length.out = dim(arr)[2])
    ),
    as.character(seq(from=ifelse(file_shape$direction["z"] > 0, file_shape$starts["z"], file_shape$ends["z"]),
                     to=ifelse(file_shape$direction["z"] > 0, file_shape$ends["z"], file_shape$starts["z"]),
                     length.out = dim(arr)[3])
    )
  )

  # Get slices
  anatomy_df <- seq_along(slices_df$slice_index) %>%
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

                  # Wrangle data to long form
                  out <- slice_array %>%
                    as_tibble %>%
                    mutate(y=rownames(slice_array)) %>%
                    gather(key = "z", value="intensity", -one_of("y")) %>%
                    mutate(x=w,
                           y=as.numeric(y),
                           z=as.numeric(z),
                           slice_axis=slice_axis,
                           slice_voxel=v,
                           slice_world=w,
                           slice_index=index) %>%
                    select(x, y, z, intensity, slice_axis, slice_voxel, slice_world, slice_index)
                },
                y={

                  # Get slice array
                  slice_array <- arr[,(v + 1),]

                  # Wrangle data to long form
                  out <- slice_array %>%
                    as_tibble %>%
                    mutate(x=rownames(slice_array)) %>%
                    gather(key = "z", value="intensity", -one_of("x")) %>%
                    mutate(x=as.numeric(x),
                           y=w,
                           z=as.numeric(z),
                           slice_axis=slice_axis,
                           slice_voxel=v,
                           slice_world=w,
                           slice_index=index) %>%
                    select(x, y, z, intensity, slice_axis, slice_voxel, slice_world, slice_index)
                },
                z={

                  # Get slice array
                  slice_array <- arr[,,(v + 1)]

                  # Wrangle data to long form
                  out <- slice_array %>%
                    as_tibble %>%
                    mutate(x=rownames(slice_array)) %>%
                    gather(key = "y", value="intensity", -one_of("x")) %>%
                    mutate(x=as.numeric(x),
                           y=as.numeric(y),
                           z=w,
                           slice_axis=slice_axis,
                           slice_voxel=v,
                           slice_world=w,
                           slice_index=index) %>%
                    select(x, y, z, intensity, slice_axis, slice_voxel, slice_world, slice_index)
                }
        )



        # Return
        return(out)
      }
    )

  # Return
  out <- list(file=file,
              anatomy_df=anatomy_df,
              slice_axis=slice_axis,
              slice_axis_coordinates=slice_axis_coordinates)
}

prepare_masked_anatomy <- function(file, mask_file, slice_axis, slice_axis_coordinates) {

  # Read data
  anatomy_data <- prepare_anatomy(file, slice_axis = slice_axis, slice_axis_coordinates = slice_axis_coordinates)
  mask_data <- prepare_anatomy(mask_file, slice_axis = slice_axis, slice_axis_coordinates = slice_axis_coordinates)

  # Merge mask data
  anatomy_and_mask_df <- anatomy_data$anatomy_df %>%
    left_join(mask_data$anatomy_df %>%
                select(x, y, z, mask_value=intensity),
              by=c("x", "y", "z")) %>%
    select(x, y, z, intensity, mask_value, everything())

  # Output
  out <- list(file=file,
              anatomy_df=anatomy_and_mask_df,
              slice_axis=slice_axis,
              slice_axis_coordinates=slice_axis_coordinates)

  return(out)

}

#############################

prepare_contours <- function(file, slice_axis, slice_axis_coordinates, contour_levels, mask_file=NULL, filter_size=NULL) {

  # Get file shape
  file_shape <- get_file_shape(file)

  # Get voxel indices
  slices_df <- get_slice_mappings(file = file, slice_axis = slice_axis, slice_axis_coordinates = slice_axis_coordinates)

  # Read in 3D file
  arr <- mincArray(mincGetVolume(file))
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

                  # Use grDevices::contourLines() to get the contours
                  # Returns a list of contour paths
                  slice_contours <- contourLines(x = seq(from=file_shape$starts["y"],
                                                         to=file_shape$ends["y"],
                                                         length.out = dim(arr)[[2]]),
                                                 y = seq(from=file_shape$starts["z"],
                                                         to=file_shape$ends["z"],
                                                         length.out = dim(arr)[[3]]),
                                                 z = slice_array,
                                                 levels = contour_levels)

                  # Wrangle data to long form
                  out <- 1:length(slice_contours) %>%
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
                           slice_index=index) %>%
                    select(level, obj, x, y, z, slice_axis, slice_voxel, slice_world, slice_index)
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

                  # Use grDevices::contourLines() to get the contours
                  # Returns a list of contour paths
                  slice_contours <- contourLines(x = seq(from=file_shape$starts["x"],
                                                         to=file_shape$ends["x"],
                                                         length.out = dim(arr)[[1]]),
                                                 y = seq(from=file_shape$starts["z"],
                                                         to=file_shape$ends["z"],
                                                         length.out = dim(arr)[[3]]),
                                                 z = slice_array,
                                                 levels = contour_levels)

                  # Wrangle data to long form
                  out <- 1:length(slice_contours) %>%
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
                           slice_index=index) %>%
                    select(level, obj, x, y, z, slice_axis, slice_voxel, slice_world, slice_index)
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

                  # Use grDevices::contourLines() to get the contours
                  # Returns a list of contour paths
                  slice_contours <- contourLines(x = seq(from=file_shape$starts["x"],
                                                         to=file_shape$ends["x"],
                                                         length.out = dim(arr)[[1]]),
                                                 y = seq(from=file_shape$starts["y"],
                                                         to=file_shape$ends["y"],
                                                         length.out = dim(arr)[[2]]),
                                                 z = slice_array,
                                                 levels = contour_levels)

                  # Wrangle data to long form
                  out <- 1:length(slice_contours) %>%
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
                           slice_index=index) %>%
                    select(level, obj, x, y, z, slice_axis, slice_voxel, slice_world, slice_index)
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
      filter(length(level) >= filter_size) %>%
      ungroup()
  }

  # Return
  out <- list(file=file,
              contours_df=contours_df,
              slice_axis=slice_axis,
              slice_axis_coordinates=slice_axis_coordinates,
              contour_levels=contour_levels,
              mask_file=mask_file,
              filter_size=filter_size)
}

#############################

# Visualize grid, and underlying anatomy if df provided
visualize_grid <- function(grid_object,
                           flip_axes=F,
                           filter_grid="none",
                           interpolate_anatomy=F,
                           grid_color='#FF6666', grid_size=0.25, grid_alpha=0.8,
                           box_color='orange', box_size=1, box_alpha=1,
                           point_color='red', point_size=0.05, point_alpha=0.2,
                           anatomy_low="auto", anatomy_high="auto") {

  anatomy_df <- prepare_anatomy(file = grid_object$file,
                                slice_axis=grid_object$slice_axis,
                                slice_axis_coordinates=grid_object$slice_axis_coordinate)$anatomy_df

  # Get axes
  all_axes <- c("x", "y", "z")
  grid_and_line_axes <- sort(setdiff(all_axes, grid_object$slice_axis))

  # Flip axes if required
  if (flip_axes) {
    horizontal_axis <- grid_and_line_axes[[2]]
    vertical_axis <- grid_and_line_axes[[1]]
  } else {
    horizontal_axis <- grid_and_line_axes[[1]]
    vertical_axis <- grid_and_line_axes[[2]]
  }

  # Filter grid if required
  slice_axis <- grid_object$slice_axis
  if (filter_grid=="positive") {
    fg_grid <- grid_object$grid %>% filter(!!sym(slice_axis) >= grid_object$slice_axis_coordinate)
    fg_box <- grid_object$box %>% filter(!!sym(slice_axis) >= grid_object$slice_axis_coordinate)
  } else if (filter_grid=="negative") {
    fg_grid <- grid_object$grid %>% filter(!!sym(slice_axis) <= grid_object$slice_axis_coordinate)
    fg_box <- grid_object$box %>% filter(!!sym(slice_axis) <= grid_object$slice_axis_coordinate)
  } else {
    fg_grid <- grid_object$grid
    fg_box <- grid_object$box
  }

  # Get anatomy scale
  if (!is.numeric(anatomy_low)) {
    anatomy_low <- quantile(anatomy_df$intensity, c(0.50))
  }
  if (!is.numeric(anatomy_high)) {
    anatomy_high <- quantile(anatomy_df$intensity, c(0.99))
  }

  print(
    plt <- fg_grid %>%
      ggplot(aes_string(x=horizontal_axis, y=vertical_axis)) +
      geom_raster(aes(fill=intensity), data=anatomy_df, interpolate = interpolate_anatomy) +
      geom_point(size=point_size, color=point_color, alpha=point_alpha) +
      geom_path(aes(group=interaction(grid_axis, gridline)), color=grid_color, size=grid_size, alpha=grid_alpha) +
      geom_path(aes(group=interaction(grid_axis, gridline)), color=box_color, size=box_size, alpha=box_alpha, data=fg_box) +
      coord_fixed(ratio=1) +
      scale_fill_gradient(low="black", high = "white", limits=c(anatomy_low, anatomy_high), oob=squish) +
      xlab(glue("{horizontal_axis}-coordinate (mm)")) +
      ylab(glue("{vertical_axis}-coordinate (mm)")) +
      labs(title=glue("{horizontal_axis}{vertical_axis}-slice at {grid_object$slice_axis}-coordinate = {round(grid_object$slice_axis_coordinate,3)} mm")) +
      theme_bw()
  )

  return(plt)
}

#############################

write_tag <- function(df, outfile, clobber=F) {

  # Check file exists
  if (file.exists(outfile)) {
    print(glue("Tag file exists: {outfile}"))
    if (clobber) {
      print(glue("Overwriting tag file..."))
    } else {
      stop("Set clobber=TRUE to overwrite existing tag file")
    }
  } else {
    print(glue("Creating tag file: {outfile}"))
  }

  # Open connection to file
  fcon <- file(outfile, open = "wt")

  # Header
  writeLines("MNI Tag Point File", con = fcon, sep = "\n")
  writeLines("Volumes = 1;", con = fcon, sep = "\n")
  writeLines("", con = fcon, sep = "\n")
  writeLines("Points =", con = fcon, sep = "\n")

  # Write points
  prog <- txtProgressBar(max=nrow(df), style=3)
  for (i in 1:nrow(df)) {
    writeLines(paste0(" ", df$x[i], " ", df$y[i], " ", df$z[i], " ", "\"\""), con = fcon, sep = "\n")
    setTxtProgressBar(prog, i)
  }
  writeLines(";", con = fcon, sep = "\n")

  # Close file connection
  close(fcon)

  # Done
  print(glue("Done writing tag file: {outfile}"))
}

#############################

read_tag <- function(tagfile, transformed=F) {

  n_lines <- as.integer(strsplit(system(glue("wc -l {tagfile}"), intern = T), " ")[[1]][1])
  if (transformed) {
    n_points <- n_lines - 4
    skip_lines <- 4
    tag_col_names <- c("V1", "x", "y", "z", "V2", "V3", "V4", "V5")
  } else {
    n_points <- n_lines - 5
    skip_lines <- 4
    tag_col_names <- c("V1", "x", "y", "z", "V2")
  }


  out <- read.table(tagfile,
                    header=F,
                    col.names = tag_col_names,
                    sep=" ",
                    nrows = n_points,
                    skip=skip_lines) %>%
    select(x,y,z)

  return(out)
}

transform_tag <- function(tagfile, xfmfile, outfile, invert=F) {
  cmd <- glue("transform_tags {tagfile} {xfmfile} {outfile}")
  if (invert) {
    cmd <- glue("{cmd} invert")
  }
  system(cmd)
}

transform_grids <- function(base_grid,
                            xfmfile,
                            invert=F,
                            transformed_space=NA,
                            transformed_file=NA,
                            tmpdir="/tmp") {

  # Tempfiles
  prexfm_grid_tagfile <- tempfile(pattern = "transform_grids-grid-pre-", tmpdir = tmpdir, fileext = ".tag")
  prexfm_box_tagfile <- tempfile(pattern = "transform_grids-box-pre-", tmpdir = tmpdir, fileext = ".tag")
  postxfm_grid_tagfile <- tempfile(pattern = "transform_grids-grid-post-", tmpdir = tmpdir, fileext = ".tag")
  postxfm_box_tagfile <- tempfile(pattern = "transform_grids-box-post-", tmpdir = tmpdir, fileext = ".tag")

  # Write out tags
  write_tag(base_grid$grid, prexfm_grid_tagfile, clobber = T)
  write_tag(base_grid$box, prexfm_box_tagfile, clobber = T)

  # Transform tags
  transform_tag(tagfile = prexfm_grid_tagfile, xfmfile = xfmfile, outfile = postxfm_grid_tagfile, invert = invert)
  transform_tag(tagfile = prexfm_box_tagfile, xfmfile = xfmfile, outfile = postxfm_box_tagfile, invert = invert)

  # Read tags
  transformed_grid_tags <- read_tag(postxfm_grid_tagfile, transformed = T)
  transformed_box_tags <- read_tag(postxfm_box_tagfile, transformed = T)

  # Cleanup tags
  file.remove(c(prexfm_grid_tagfile, prexfm_box_tagfile, postxfm_grid_tagfile, postxfm_box_tagfile))

  # Label tags
  transformed_grid_df <- base_grid$grid %>%
    mutate(grid_space=transformed_space) %>%
    select(-x, -y, -z) %>%
    cbind(transformed_grid_tags)

  transformed_box_df <- base_grid$box %>%
    mutate(grid_space=transformed_space) %>%
    select(-x, -y, -z) %>%
    cbind(transformed_box_tags)

  # Output
  transformation_info=list(
    file=transformed_file,
    xfmfile=xfmfile,
    invert=invert,
    pre_xfm_file=base_grid$file,
    pre_xfm_grid_space=base_grid$grid_space,
    pre_xfm_slice_axis=base_grid$slice_axis,
    pre_xfm_slice_axis_coordinate=base_grid$slice_axis_coordinate
  )

  out <- list(grid=transformed_grid_df,
              box=transformed_box_df,
              file=transformed_file,
              grid_space=transformed_space,
              slice_axis=base_grid$slice_axis,
              slice_axis_coordinate=mean(transformed_grid_df[[base_grid$slice_axis]]),
              transformation_info=transformation_info)

  # Return
  return(out)

}


# Input: nlin coronal slice, zoomed coordinates on nlin slice
# Output:
get_pipeline_grids <- function(slice_axis, nlin_slice_axis_coordinate,
                               native_file, lsq6_file, nlin_file,
                               study_template_file, study_mask_file,
                               native_to_lsq6_xfm, lsq6_to_nlin_xfm,
                               nlin_abs_jd_file,
                               nlin_contour_levels,
                               grid_padding = 0.5, line_points = 200, grid_spacing = 0.2,
                               highres_grid_lines=12,
                               tmpdir="/tmp", plot_progress=T) {

  #############
  # Data preparation
  #############

  # Get axes
  cat(glue("[{Sys.time()}] Getting axis information \n", .trim=F))
  all_axes <- c("x", "y", "z")
  grid_and_line_axes <- sort(setdiff(all_axes, slice_axis))

  # Get nlin grid and transform to lsq6 to get coordinates there
  cat(glue("[{Sys.time()}] Determining optimal lsq6 slice coordinate\n", .trim=F))

  cat(glue("[{Sys.time()}] + Creating base grid on nlin slice \n", .trim=F))
  nlin_gs <- get_grid_sequence(file = nlin_file, grid_padding = 0, line_points = 50, grid_spacing = 0.2)
  nlin_grid <- get_base_grid(nlin_gs, grid_space = "nlin", slice_axis = slice_axis, slice_axis_coordinate = nlin_slice_axis_coordinate)
  if (plot_progress) {
    visualize_grid(nlin_grid, flip_axes = F)
  }

  # Transform to lsq6
  cat(glue("[{Sys.time()}] + Transforming nlin grid on lsq6 \n", .trim=F))
  nlin_grid_on_lsq6 <- transform_grids(nlin_grid,
                                       xfmfile = lsq6_to_nlin_xfm,
                                       invert = T,
                                       transformed_space = "lsq6",
                                       transformed_file = lsq6_file,
                                       tmpdir = tmpdir)
  if (plot_progress) {
    visualize_grid(nlin_grid_on_lsq6, filter_grid = "none", interpolate_anatomy = F)
  }

  cat(glue("[{Sys.time()}] + Getting lsq6 slice coordinate \n", .trim=F))
  lsq6_slice_axis_coordinate <- nlin_grid_on_lsq6$slice_axis_coordinate

  #############
  # Get study template data
  #############

  cat(glue("[{Sys.time()}] Getting study data \n", .trim=F))

  cat(glue("[{Sys.time()}] + Getting study template data \n", .trim=F))
  study_template_data <- prepare_masked_anatomy(study_template_file,
                                                mask_file = study_mask_file,
                                                slice_axis = slice_axis,
                                                slice_axis_coordinates = nlin_slice_axis_coordinate)

  cat(glue("[{Sys.time()}] + Getting study mask data \n", .trim=F))
  study_mask_data <- prepare_anatomy(study_mask_file,
                                     slice_axis = slice_axis,
                                     slice_axis_coordinates = nlin_slice_axis_coordinate)

  cat(glue("[{Sys.time()}] + Getting full contour set \n", .trim=F))
  study_full_contour_data <- prepare_contours(study_template_file,
                                              slice_axis = slice_axis,
                                              slice_axis_coordinates = nlin_slice_axis_coordinate,
                                              contour_levels = nlin_contour_levels)

  cat(glue("[{Sys.time()}] + Getting masked contour set \n", .trim=F))
  study_masked_contour_data <- prepare_contours(study_template_file,
                                                mask_file = study_mask_file,
                                                slice_axis = slice_axis,
                                                slice_axis_coordinates = nlin_slice_axis_coordinate,
                                                contour_levels = nlin_contour_levels)

  #############
  # Get Jacobian determinant data
  #############

  # Get JD data
  cat(glue("[{Sys.time()}] Getting Jacobian determinant data \n", .trim=F))
  jd_data <- prepare_masked_anatomy(nlin_abs_jd_file,
                                    mask_file = study_mask_file,
                                    slice_axis = slice_axis,
                                    slice_axis_coordinates = nlin_slice_axis_coordinate)

  # Get x,y,z coordinates where JD is highest within the brain
  cat(glue("[{Sys.time()}] + Getting JD data within brain \n", .trim=F))
  jd_intensities <- jd_data$anatomy_df %>%
    left_join(study_mask_data$anatomy_df %>%
                select(x, y, z, is_brain=intensity),
              by=c("x", "y", "z")) %>%
    filter(is_brain==1) %>%
    arrange(desc(intensity)) %>%
    select(x, y, z, logjd=intensity)

  cat(glue("[{Sys.time()}] + Getting JD extremes \n", .trim=F))
  jd_expansion_coordinates <- jd_intensities %>%
    head(1) %>%
    unlist()
  jd_contraction_coordinates <- jd_intensities %>%
    tail(1) %>%
    unlist()

  cat(glue("[{Sys.time()}] + Merging into dataframe \n", .trim=F))
  jd_coordinates_nlin <- rbind(jd_intensities %>%
                                 head(1) %>%
                                 mutate(grid_space="nlin",
                                        note="expansion"),
                               jd_intensities %>%
                                 tail(1) %>%
                                 mutate(grid_space="nlin",
                                        note="contraction")
  )


  #############
  # Get anatomy and grid data for lsq6
  #############

  cat(glue("[{Sys.time()}] Getting lsq6 data \n", .trim=F))

  cat(glue("[{Sys.time()}] + Getting lsq6 grid \n", .trim=F))
  lsq6_gs <- get_grid_sequence(file = lsq6_file, grid_padding = grid_padding, line_points = line_points, grid_spacing = grid_spacing)
  lsq6_grid <- get_base_grid(lsq6_gs, grid_space = "lsq6", slice_axis = slice_axis, slice_axis_coordinate = lsq6_slice_axis_coordinate)

  cat(glue("[{Sys.time()}] + Getting lsq6 anatomy data \n", .trim=F))
  lsq6_data <- prepare_anatomy(file = lsq6_grid$file,
                               slice_axis=lsq6_grid$slice_axis,
                               slice_axis_coordinates=lsq6_grid$slice_axis_coordinate)
  if (plot_progress) {
    visualize_grid(lsq6_grid, flip_axes = F)
  }

  #############
  # Get anatomy and grid data for native
  #############

  cat(glue("[{Sys.time()}] Getting native data \n", .trim=F))

  cat(glue("[{Sys.time()}] + Getting lsq6 grid on native \n", .trim=F))
  lsq6_grid_on_native <- transform_grids(lsq6_grid,
                                         xfmfile = native_to_lsq6_xfm,
                                         invert = T,
                                         transformed_space = "native",
                                         transformed_file = native_file,
                                         tmpdir = tmpdir)

  cat(glue("[{Sys.time()}] + Getting native anatomy data \n", .trim=F))
  native_data <- prepare_anatomy(file = lsq6_grid_on_native$file,
                                 slice_axis=lsq6_grid_on_native$slice_axis,
                                 slice_axis_coordinates=lsq6_grid_on_native$slice_axis_coordinate)
  if (plot_progress) {
    visualize_grid(lsq6_grid_on_native, filter_grid = "positive")
  }

  #############
  # Get anatomy and grid data for nlin
  #############

  cat(glue("[{Sys.time()}] Getting nlin data \n", .trim=F))

  cat(glue("[{Sys.time()}] + Getting lsq6 grid on nlin \n", .trim=F))
  lsq6_grid_on_nlin <- transform_grids(lsq6_grid,
                                       xfmfile = lsq6_to_nlin_xfm,
                                       invert = F,
                                       transformed_space = "nlin",
                                       transformed_file = nlin_file,
                                       tmpdir = tmpdir)

  cat(glue("[{Sys.time()}] + Getting nlin anatomy data \n", .trim=F))
  nlin_data <- prepare_anatomy(file = lsq6_grid_on_nlin$file,
                               slice_axis=lsq6_grid_on_nlin$slice_axis,
                               slice_axis_coordinates=lsq6_grid_on_nlin$slice_axis_coordinate)

  if (plot_progress) {
    visualize_grid(lsq6_grid_on_nlin, filter_grid = "none")
  }

  #############
  # Get JD point mapping in lsq6 space
  #############

  cat(glue("[{Sys.time()}] Mapping JD extrema to lsq6 space \n", .trim=F))
  cat(glue("[{Sys.time()}] + Mapping via tags \n", .trim=F))

  # Transform points to lsq6
  prexfm_point_tagfile <- tempfile(pattern = "transform_grids-point-pre-", tmpdir = tmpdir, fileext = ".tag")
  postxfm_point_tagfile <- tempfile(pattern = "transform_grids-point-post-", tmpdir = tmpdir, fileext = ".tag")

  # Write out tags
  write_tag(jd_coordinates_nlin, prexfm_point_tagfile, clobber = T)

  # Transform tags
  transform_tag(tagfile=prexfm_point_tagfile, xfmfile = lsq6_to_nlin_xfm, outfile = postxfm_point_tagfile, invert = T)

  # Read tags
  transformed_point_tags <- read_tag(postxfm_point_tagfile, transformed = T)

  # Cleanup tags
  file.remove(c(prexfm_point_tagfile, postxfm_point_tagfile))

  # Label tags
  cat(glue("[{Sys.time()}] + Merging transformed tags with row information \n", .trim=F))
  jd_coordinates_lsq6 <- jd_coordinates_nlin %>%
    mutate(grid_space="lsq6") %>%
    select(-x, -y, -z) %>%
    cbind(transformed_point_tags)

  #############
  # Get voxel resolution grid sequences for lsq6
  #############

  cat(glue("[{Sys.time()}] Getting high resolution grid sequence in lsq6 space \n", .trim=F))

  lsq6_gs_highres_expansion <- get_grid_sequence_around_point(file = lsq6_file,
                                                              point = jd_coordinates_lsq6[jd_coordinates_lsq6$note=="expansion",c("x", "y", "z")],
                                                              grid_spacing = mean(get_file_shape(nlin_file)$steps[grid_and_line_axes]),
                                                              grid_lines = highres_grid_lines,
                                                              line_points = highres_grid_lines)

  lsq6_gs_highres_contraction <- get_grid_sequence_around_point(file = lsq6_file,
                                                                point = jd_coordinates_lsq6[jd_coordinates_lsq6$note=="contraction",c("x", "y", "z")],
                                                                grid_spacing = mean(get_file_shape(nlin_file)$steps[grid_and_line_axes]),
                                                                grid_lines = highres_grid_lines,
                                                                line_points = highres_grid_lines)



  ############
  # Get voxel resolution grid data for lsq6
  #############

  cat(glue("[{Sys.time()}] Getting high resolution grid data in lsq6 space \n", .trim=F))

  lsq6_highres_grid_expansion <- get_base_grid(lsq6_gs_highres_expansion,
                                               grid_space = "lsq6",
                                               slice_axis = slice_axis,
                                               slice_axis_coordinate = jd_coordinates_lsq6[jd_coordinates_lsq6$note=="expansion",slice_axis])
  if (plot_progress) {
    visualize_grid(lsq6_highres_grid_expansion)
  }

  lsq6_highres_grid_contraction <- get_base_grid(lsq6_gs_highres_contraction,
                                                 grid_space = "lsq6",
                                                 slice_axis = slice_axis,
                                                 slice_axis_coordinate = jd_coordinates_lsq6[jd_coordinates_lsq6$note=="contraction",slice_axis])
  if (plot_progress) {
    visualize_grid(lsq6_highres_grid_contraction)
  }



  #############
  # Get voxel resolution grid data for nlin
  #############

  cat(glue("[{Sys.time()}] Getting high resolution grid sequence in nlin space \n", .trim=F))

  # Get high res grids at expansion and contraction points in lsq6 space
  lsq6_highres_grid_expansion_on_nlin <- transform_grids(lsq6_highres_grid_expansion,
                                                         xfmfile = lsq6_to_nlin_xfm,
                                                         invert = F,
                                                         transformed_space = "nlin",
                                                         transformed_file = nlin_file,
                                                         tmpdir = tmpdir)
  if (plot_progress) {
    visualize_grid(lsq6_highres_grid_expansion_on_nlin)
  }


  lsq6_highres_grid_contraction_on_nlin <- transform_grids(lsq6_highres_grid_contraction,
                                                           xfmfile = lsq6_to_nlin_xfm,
                                                           invert = F,
                                                           transformed_space = "nlin",
                                                           transformed_file = nlin_file,
                                                           tmpdir = tmpdir)
  if (plot_progress) {
    visualize_grid(lsq6_highres_grid_contraction_on_nlin)
  }


  #############
  # Output
  #############

  cat(glue("[{Sys.time()}] Computing outputs \n", .trim=F))

  inputs <- list(
    slice_axis=slice_axis,
    nlin_slice_axis_coordinate=nlin_slice_axis_coordinate,
    native_file=native_file,
    lsq6_file=lsq6_file,
    nlin_file=nlin_file,
    study_template_file=study_template_file,
    study_mask_file=study_mask_file,
    native_to_lsq6_xfm=native_to_lsq6_xfm,
    lsq6_to_nlin_xfm=lsq6_to_nlin_xfm,
    nlin_abs_jd_file=nlin_abs_jd_file,
    grid_padding=grid_padding,
    line_points=line_points,
    grid_spacing=grid_spacing,
    highres_grid_lines=highres_grid_lines,
    tmpdir=tmpdir
  )

  lsq6_outputs=list(
    grid=lsq6_grid,
    highres_grid_expansion=lsq6_highres_grid_expansion,
    highres_grid_contraction=lsq6_highres_grid_contraction,
    anatomy=lsq6_data
  )

  native_outputs=list(
    grid=lsq6_grid_on_native,
    anatomy=native_data
  )

  nlin_outputs=list(
    grid=lsq6_grid_on_nlin,
    highres_grid_expansion=lsq6_highres_grid_expansion_on_nlin,
    highres_grid_contraction=lsq6_highres_grid_contraction_on_nlin,
    anatomy=nlin_data
  )

  study_outputs=list(
    template=study_template_data,
    mask=study_mask_data,
    contours=study_full_contour_data,
    contours_masked=study_masked_contour_data
  )

  jd_outputs=list(
    jd_data=jd_data,
    jd_coordinates_lsq6=jd_coordinates_lsq6,
    jd_coordinates_nlin=jd_coordinates_nlin
  )

  cat(glue("[{Sys.time()}] Binding outputs \n", .trim=F))

  out <- list(
    inputs=inputs,
    native=native_outputs,
    lsq6=lsq6_outputs,
    nlin=nlin_outputs,
    study=study_outputs,
    jd=jd_outputs
  )

  # Return
  return(out)

}



