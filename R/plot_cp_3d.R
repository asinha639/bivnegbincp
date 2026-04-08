#' Interactive 3D Plot of Bivariate Time Series and Change Points
#'
#' Generates an interactive 3D plot using plotly to visualize a bivariate time
#' series, with detected change points highlighted in a polished 3D display.
#' The first two axes correspond to the two observed variables, while the third
#' axis represents time.
#'
#' @param x1 A numeric vector for the first variable.
#' @param x2 A numeric vector for the second variable.
#' @param change_points A numeric vector of detected change point indices.
#' @param var1_name Character string for the first variable axis label.
#'   Default is `"Variable 1"`.
#' @param var2_name Character string for the second variable axis label.
#'   Default is `"Variable 2"`.
#' @param time_name Character string for the time axis label.
#'   Default is `"Time"`.
#' @param theme Character string specifying the visual theme. Either `"dark"`
#'   or `"light"`. Default is `"dark"`.
#' @param point_size Numeric marker size for the main data points.
#'   Default is `5`.
#' @param cp_size Numeric marker size for the highlighted change points.
#'   Default is `7`.
#'
#' @return A plotly interactive 3D plot object.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x1 <- c(rpois(20, 20), rpois(15, 55), rpois(15, 95))
#' x2 <- c(rpois(20, 10), rpois(15, 45), rpois(15, 18))
#' cp <- c(20, 35)
#'
#' plot_cp_3d(
#'   x1 = x1,
#'   x2 = x2,
#'   change_points = cp,
#'   var1_name = "Cases",
#'   var2_name = "Deaths",
#'   time_name = "Day",
#'   theme = "dark"
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @export
plot_cp_3d <- function(x1,
                       x2,
                       change_points,
                       var1_name = "Variable 1",
                       var2_name = "Variable 2",
                       time_name = "Time",
                       theme = "dark",
                       point_size = 5,
                       cp_size = 7) {
  theme <- match.arg(theme, c("dark", "light"))
  
  if (!is.numeric(x1) || !is.numeric(x2)) {
    stop("'x1' and 'x2' must be numeric vectors.")
  }
  
  if (length(x1) != length(x2)) {
    stop("'x1' and 'x2' must have the same length.")
  }
  
  if (length(x1) < 2) {
    stop("At least two observations are required.")
  }
  
  if (!is.numeric(change_points)) {
    stop("'change_points' must be a numeric vector.")
  }
  
  n <- length(x1)
  timeline <- seq_len(n)
  
  d1 <- data.frame(
    x1 = x1,
    x2 = x2,
    timeline = timeline
  )
  
  change_points <- unique(as.integer(change_points))
  change_points <- change_points[!is.na(change_points)]
  change_points <- change_points[change_points >= 1 & change_points <= n]
  
  d_cp <- data.frame(
    x = x1[change_points],
    y = x2[change_points],
    z = change_points
  )
  
  if (theme == "dark") {
    paper_bg <- "#0b1020"
    plot_bg  <- "#111827"
    grid_col <- "rgba(255,255,255,0.10)"
    zero_col <- "rgba(255,255,255,0.20)"
    font_col <- "#e5e7eb"
    point_col <- "#38bdf8"
    point_line_col <- "#e0f2fe"
    cp_col <- "#f43f5e"
    cp_line_col <- "#fb7185"
  } else {
    paper_bg <- "#f8fafc"
    plot_bg  <- "#ffffff"
    grid_col <- "rgba(15,23,42,0.12)"
    zero_col <- "rgba(15,23,42,0.25)"
    font_col <- "#0f172a"
    point_col <- "#2563eb"
    point_line_col <- "#bfdbfe"
    cp_col <- "#dc2626"
    cp_line_col <- "#ef4444"
  }
  
  p <- plotly::plot_ly()
  
  p <- p %>%
    plotly::add_trace(
      data = d1,
      x = ~x1,
      y = ~x2,
      z = ~timeline,
      type = "scatter3d",
      mode = "markers",
      name = "Series",
      hovertemplate = paste0(
        var1_name, ": %{x}<br>",
        var2_name, ": %{y}<br>",
        time_name, ": %{z}<extra></extra>"
      ),
      marker = list(
        size = point_size,
        color = point_col,
        opacity = 0.95,
        line = list(color = point_line_col, width = 1.2),
        symbol = "circle"
      )
    )
  
  if (nrow(d_cp) > 0) {
    p <- p %>%
      plotly::add_trace(
        data = d_cp,
        x = ~x,
        y = ~y,
        z = ~z,
        type = "scatter3d",
        mode = "markers",
        name = "Change points",
        hovertemplate = paste0(
          "Change point<br>",
          var1_name, ": %{x}<br>",
          var2_name, ": %{y}<br>",
          time_name, ": %{z}<extra></extra>"
        ),
        marker = list(
          size = cp_size,
          color = cp_col,
          opacity = 1,
          line = list(color = "#ffffff", width = 2),
          symbol = "circle"
        )
      ) %>%
      plotly::add_trace(
        data = d_cp,
        x = ~x,
        y = ~y,
        z = ~z,
        type = "scatter3d",
        mode = "lines",
        showlegend = FALSE,
        hoverinfo = "skip",
        line = list(
          color = cp_line_col,
          width = 5
        )
      )
  }
  
  p %>%
    plotly::layout(
      paper_bgcolor = paper_bg,
      plot_bgcolor = plot_bg,
      font = list(
        family = "Arial, sans-serif",
        size = 14,
        color = font_col
      ),
      margin = list(l = 0, r = 0, b = 0, t = 30),
      legend = list(
        orientation = "h",
        x = 0.02,
        y = 0.98,
        bgcolor = "rgba(0,0,0,0)"
      ),
      scene = list(
        bgcolor = plot_bg,
        camera = list(
          eye = list(x = 1.7, y = 1.7, z = 1.1)
        ),
        xaxis = list(
          title = list(text = var1_name),
          showbackground = TRUE,
          backgroundcolor = plot_bg,
          gridcolor = grid_col,
          zerolinecolor = zero_col,
          color = font_col
        ),
        yaxis = list(
          title = list(text = var2_name),
          showbackground = TRUE,
          backgroundcolor = plot_bg,
          gridcolor = grid_col,
          zerolinecolor = zero_col,
          color = font_col
        ),
        zaxis = list(
          title = list(text = time_name),
          showbackground = TRUE,
          backgroundcolor = plot_bg,
          gridcolor = grid_col,
          zerolinecolor = zero_col,
          color = font_col
        )
      )
    )
}