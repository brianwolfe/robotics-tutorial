
nonlinear_comparison = filter_plotting.init_nonlinear_comparison()

d3.select('#particlefilter')
    .datum([
    {
      key: 'particlefilter',
      mu: [[0], [0], [0]],
      sigma: [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 0.5]],
      type: 'particle',
    }
    ]).call(nonlinear_comparison);


d3.select('#comparison')
    .datum([
    {
      key: 'particlefilter',
      mu: [[0], [0], [0]],
      sigma: [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 0.5]],
      type: 'particle',
    },
    {
      key: 'unscentedfilter',
      mu: [[0], [0], [0]],
      sigma: [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 0.5]],
      type: 'ukf',
    }
    ]).call(nonlinear_comparison);
