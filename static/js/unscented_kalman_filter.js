

nonlinear_comparison = filter_plotting.init_nonlinear_comparison()

d3.select('#singlemove')
    .datum([
    {
      key: "initial",
      mu: [[0], [0], [0]],
      sigma: [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
    }
    ])
    .call(nonlinear_comparison);
