library(ggplot2)

a <- 0
b <- 1
N <- 5
Delta <- (b - a) / N
t_grid <- a + ((2 * seq_len(N) - 1) / (2 * N)) * (b - a)
z <- c(1.10, 1.70, 1.35, 1.95, 1.45)

embedding <- data.frame(
  t = c(a, t_grid, b),
  value = c(z[1], z, z[N])
)

observations <- data.frame(
  t = t_grid,
  value = z,
  n = seq_len(N)
)

interval_labels <- data.frame(
  x = c((a + t_grid[1]) / 2, (t_grid[1] + t_grid[2]) / 2),
  y = min(z) - 0.30,
  label = c("Delta/2", "Delta")
)

interval_arrows <- data.frame(
  x = c(a, t_grid[1]),
  xend = c(t_grid[1], t_grid[2]),
  y = min(z) - 0.22,
  yend = min(z) - 0.22
)

p <- ggplot() +
  geom_vline(
    xintercept = c(a, t_grid, b),
    colour = "grey78",
    linewidth = 0.35,
    linetype = "22"
  ) +
  geom_segment(
    data = interval_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(ends = "both", length = unit(0.06, "in")),
    linewidth = 0.35,
    colour = "grey35"
  ) +
  geom_text(
    data = interval_labels,
    aes(x = x, y = y, label = label),
    parse = TRUE,
    size = 3.4,
    vjust = 1.1,
    colour = "grey25"
  ) +
  geom_line(
    data = embedding,
    aes(x = t, y = value),
    colour = "#D97706",
    linewidth = 1.1
  ) +
  geom_point(
    data = observations,
    aes(x = t, y = value),
    colour = "#2563EB",
    fill = "white",
    shape = 21,
    stroke = 1.0,
    size = 3.2
  ) +
  geom_text(
    data = observations,
    aes(x = t, y = value, label = paste0("z[", n, "]")),
    parse = TRUE,
    nudge_y = 0.13,
    size = 3.6,
    colour = "#1E3A8A"
  ) +
  annotate(
    "text",
    x = 0.82,
    y = 1.86,
    label = "\u2130[z](t)",
    parse = TRUE,
    family = "STIX Two Math",
    size = 5.2,
    colour = "#92400E"
  ) +
  scale_x_continuous(
    breaks = c(a, t_grid, b),
    labels = parse(text = c("a", paste0("t[", seq_len(N), "]"), "b")),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = NULL,
    limits = c(min(z) - 0.42, max(z) + 0.28)
  ) +
  labs(
    x = NULL,
    title = "Piecewise linear embedding",
    subtitle = "Midpoint grid observations mapped to a continuous function"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, colour = "grey30"),
    panel.grid = element_blank(),
    axis.text.x = element_text(colour = "grey20"),
    plot.margin = margin(8, 12, 8, 12)
  )

ggsave(
  filename = "figures/linear_embedding_figure.pdf",
  plot = p,
  width = 6.4,
  height = 3.3,
  device = cairo_pdf
)

