#' @importFrom magrittr "%>%"

clip <- function(m, a, b) {
  ifelse(m <= a, a, ifelse(m >= b, b, m))
}

compute_tetrahedron_volume <- function(face, point) {
  n_compare <- face[[1]]
  n <- pracma::cross(face[[3]] - face[[2]], face[[4]] - face[[2]])
  # Qhull implementation does not create faces in counter clockwise
  # order so we might need to correct the sign of the face normals.
  if (sign(n[[1]]) != sign(n_compare[[1]])) {
    n <- -n
  }

  abs(pracma::dot(n, point - face[[4]])) / 6.0
}

# remove_one_edge_by_finding_smallest_adding_volume_with_test_condition
remove_one_edge <- function(faces, normals, edges, vertices, v2f) {
  temp_list1 <- list()
  temp_list2 <- list()
  count <- 0

  for (i in seq(nrow(edges))) {
    vertex1 <- edges[i, "i1"]
    vertex2 <- edges[i, "i2"]
    id1 <- edges[i, "alt_i1"]
    id2 <- edges[i, "alt_i2"]

    related_face_ids <- v2f[
      v2f[, "p"] == vertex1 | v2f[, "p"] == vertex2,
      "face_id"
    ] %>%
      unique

    old_face_list <- list()

    c <- c(0, 0, 0)
    a <- list()
    b <- list()

    for (rel_face_id in related_face_ids) {
      face <- faces[rel_face_id, c("p1", "p2", "p3")]
      p1 <- vertices[toString(face[["p1"]]), c("r", "g", "b")]
      p2 <- vertices[toString(face[["p2"]]), c("r", "g", "b")]
      p3 <- vertices[toString(face[["p3"]]), c("r", "g", "b")]
      n <- normals[rel_face_id, c("r", "g", "b")] %>% unname()
      old_face_list[[length(old_face_list) + 1]] <- list(n, p1, p2, p3)

      a[[length(a) + 1]] <- n
      b[[length(b) + 1]] <- pracma::dot(n, p1)
      c <- c + n
    }

    a <- -matrix(unlist(a), ncol = 3, byrow = TRUE)
    b <- -unlist(b)

    res <- Rglpk::Rglpk_solve_LP(
      c, a, dir = rep("<=", length(b)), b,
      bounds = list(
        lower = list(ind = 1:3, val = rep(-Inf, 3)),
        upper = list(ind = 1:3, val = rep(Inf, 3))
      ),
      control = list(
        # "presolve" = TRUE,
        # "verbose" = TRUE,
        # "canonicalize_status" = FALSE
      )
    )
    if (res$status == 0) {  # optimal
      newpoint <- res$solution %>% drop()
      volume <- 0
      for (each_face in old_face_list) {
        volume <- volume + compute_tetrahedron_volume(each_face, newpoint)
      }
      temp_list1[[length(temp_list1) + 1]] <- c(
        count, volume, vertex1, vertex2, id1, id2
      )
      temp_list2[[length(temp_list2) + 1]] <- newpoint
      count <- count + 1
    } else {
      # print("NOT optimal solution!!!")
    }
  }

  if (length(temp_list1) == 0) {
    print("length 0")
  } else {
    # Min by volume.
    ix <- which.min(sapply(temp_list1, `[[`, 2))
    min_tuple <- temp_list1[[ix]]
    final_point <- temp_list2[[ix]]

    # Vertex indices.
    v1_ind <- min_tuple[[3]]
    v2_ind <- min_tuple[[4]]

    related_faces_vertex_ind <- v2f[
      v2f[, "p"] == v1_ind | v2f[, "p"] == v2_ind,
      "face_id"
    ] %>%
      unique

    # Remove faces.
    new_faces <- faces[-related_faces_vertex_ind, ]
    # Remove dangling vertices.
    if (!is.matrix(new_faces)) {
      new_face_vertices <- new_faces[c("p1", "p2", "p3")]
    } else {
      new_face_vertices <- new_faces[, c("p1", "p2", "p3")] %>%
        t %>%
        as.vector %>%
        unique()
    }

    dangling_vertices <- setdiff(
      as.integer(rownames(vertices)),
      new_face_vertices
    )

    vertices <- vertices[setdiff(rownames(vertices), dangling_vertices), ]
    vertices <- rbind(
      vertices,
      c(
        alt_id = max(vertices[, "alt_ix"]) + 1,
        r = final_point[[1]],
        g = final_point[[2]],
        b = final_point[[3]]
      )
    )
  }

  return(vertices)
}

faces_from_hull <- function(hull) {
  faces <- cbind(seq_len(nrow(hull$hull)), hull$hull)
  colnames(faces) <- c("face_id", "p1", "p2", "p3")
  return(faces)

  faces <- tibble::tibble(
    face_id = seq_along(hull$hull[, 1]),
    p0 = hull$hull[, 1],
    p1 = hull$hull[, 2],
    p2 = hull$hull[, 3]
  ) %>%
    tidyr::pivot_longer(p0:p2, names_to = "point", values_to = "vertex_id")
  faces
}

normals_from_hull <- function(hull) {
  normals <- cbind(seq_len(nrow(hull$normals)), hull$normals)
  colnames(normals) <- c("face_id", "r", "g", "b", "w")
  return(normals)

  normals <- tibble::tibble(
    face_id = seq_along(hull$hull[, 1]),
    n0 = hull$normals[, 1],
    n1 = hull$normals[, 2],
    n2 = hull$normals[, 3]
  )
  normals
}

edges_from_hull <- function(hull, vertices) {
  edges <- list()
  cnt <- 1
  for (i in seq(nrow(hull$hull))) {
    face <- hull$hull[i, ]
    pix <- sort(c(face[[1]], face[[2]], face[[3]]))
    p1_ix <- pix[[1]]
    p2_ix <- pix[[2]]
    p3_ix <- pix[[3]]

    alt1_ix <- vertices[toString(p1_ix), "alt_ix"]
    alt2_ix <- vertices[toString(p2_ix), "alt_ix"]
    alt3_ix <- vertices[toString(p3_ix), "alt_ix"]

    edges[[cnt + 0]] <- c(p1_ix, p2_ix, alt1_ix, alt2_ix)
    edges[[cnt + 1]] <- c(p1_ix, p3_ix, alt1_ix, alt3_ix)
    edges[[cnt + 2]] <- c(p2_ix, p3_ix, alt2_ix, alt3_ix)

    cnt <- cnt + 3
  }
  edges <- edges %>% unique %>% unlist %>% matrix(ncol = 4, byrow = TRUE)
  colnames(edges) <- c("i1", "i2", "alt_i1", "alt_i2")
  return(edges)
}

vertices_from_hull <- function(hull) {
  vertex_ids <- matrix(hull$hull) %>%
    list() %>%
    unlist() %>%
    unique() %>%
    sort()

  vertices <- cbind(0:(length(vertex_ids) - 1), hull$p[vertex_ids, ])
  colnames(vertices) <- c("alt_ix", "r", "g", "b")
  rownames(vertices) <- vertex_ids
  return(vertices)

  vertices <- tibble(
    vertex_id = vertex_ids,
    r = hull$p[vertex_ids, 1],
    g = hull$p[vertex_ids, 2],
    b = hull$p[vertex_ids, 3]
  ) %>%
    rowid_to_column("alt_id") %>%
    mutate(alt_id = alt_id - 1)
  vertices
}

vertex_to_face <- function(faces) {
  lookup <- faces[, c("p1", "p2", "p3")] %>%
    t %>%
    as.vector %>%
    matrix(ncol = 1)
  lookup <- cbind(lookup, rep(faces[, "face_id"], each = 3))
  colnames(lookup) <- c("p", "face_id")
  lookup
}

# hull_simplification_determined_version
hull_simplification <- function(rgb, n_colors) {
  stop_early_rgb <- NA
  hull <- tryCatch(
    geometry::convhulln(rgb, output.options = TRUE),
    error = function(cond) {
      if (stringr::str_detect(conditionMessage(cond), "QH(6154|6013)")) {
        print("Image is essentially 1D and has no convex hull.")
        rgbnorm <- cbind(rgb, apply(rgb, 1, function(x) sqrt(sum(x^2))))
        colnames(rgbnorm) <- c("r", "g", "b", "norm")
        sorted_rgb <- rgbnorm[
          order(rgbnorm[, c("norm")],
            decreasing = TRUE
          ), c("r", "g", "b")
        ]
        stop_early_rgb <<- rbind(
          sorted_rgb[1, ],
          sorted_rgb[nrow(sorted_rgb), ]
        )
      } else {
        stop(cond)
      }
    }
  )
  if (!is.na(stop_early_rgb[[1]])) {
    colnames(stop_early_rgb) <- c("r", "g", "b")
    rownames(stop_early_rgb) <- seq_len(nrow(stop_early_rgb))
    return(stop_early_rgb)
  }

  max_loop <- 5000
  for (i in seq(max_loop)) {
    if (i %% 100 == 0) {
      print(paste("loop i:", i))
    }

    faces <- faces_from_hull(hull)
    normals <- normals_from_hull(hull)
    vertices <- vertices_from_hull(hull)
    edges <- edges_from_hull(hull, vertices)
    v2f <- vertex_to_face(faces)
    old_num <- nrow(vertices)

    v <- remove_one_edge(faces, normals, edges, vertices, v2f)
    hull <- geometry::convhulln(v[, c("r", "g", "b")], output.options = TRUE)

    if (nrow(hull$p) == old_num || nrow(hull$p) <= n_colors) {
      return(clip(hull$p, 0, 1))
    }
  }
  error(paste("Did not find", n_colors, "colors"))
}

#' Create a theme for an image
#'
#' @param path Image path
#' @param n_colors Size of theme palette
#' @return Image color theme
#' @export
get_theme_colors_from_image <- function(path, n_colors) {
  im_resized <- magick::image_read(path) %>%
    magick::image_resize("10%", filter = "Point")
  im <- im_resized %>%
    magick::image_data("rgb") %>%
    as.numeric()
  d <- dim(im)
  rgb <- reticulate::array_reshape(im, c(d[[1]] * d[[2]], 3))

  colors_rgb <- hull_simplification(rgb, n_colors = n_colors)
  colors_lms <- colorscience::RGB2LMS(colors_rgb)
  sorted_rgb <- colors_rgb[
    order(
      colors_lms[, 1],
      colors_lms[, 2],
      colors_lms[, 3],
      decreasing = TRUE
    ),
  ]
  sorted_rgb
}

#' Create themes for all images in a directory
#'
#' @param path Directory to search for images
#' @param n_colors Size of theme palette
#' @return List of image color themes
#' @export
get_theme_colors_from_dir <- function(path, n_colors, pattern = "*.jpg") {
  images <- list.files(path, pattern = pattern)
  map(
    images,
    function(x) {
      get_theme_colors_from_image(paste0(path, x), n_colors = n_colors)
    }
  )
}
