#custom asinh transform for ggplot
asinh_trans <- trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh
)