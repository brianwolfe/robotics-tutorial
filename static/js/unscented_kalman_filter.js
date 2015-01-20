

nonlinear_comparison = filter_plotting.init_nonlinear_comparison()

d3.select('#singlemove')
    .datum([
    {
      key: 'initial',
      mu: [[0], [0], [0]],
      sigma: [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 0.5]],
      type: 'ukf',
    }
    ])
    .call(nonlinear_comparison);

d3.select('#comparison')
    .datum([
    {
      key: 'UKF',
      mu: [[0], [0], [0]],
      sigma: [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 0.5]],
      type: 'ukf'
    },
    {
      key: 'particlefilter',
      mu: [[0], [0], [0]],
      sigma: [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 0.5]],
      type: 'particle',
    }
    ]).call(nonlinear_comparison);
