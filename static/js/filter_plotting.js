
filter_plotting = function() {
  var fp = {
    version: "0.0.1",
    colorwheel: ["#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"]
  };

  var range = function(start, end, step) {
    var range = [];
    var typeofStart = typeof start;
    var typeofEnd = typeof end;

    if (step === 0) {
      throw TypeError("Step cannot be zero.");
    }

    if (typeofStart == "undefined" || typeofEnd == "undefined") {
      throw TypeError("Must pass start and end arguments.");
    } else if (typeofStart != typeofEnd) {
      throw TypeError("Start and end arguments must be of same type.");
    }

    typeof step == "undefined" && (step = 1);

    if (end < start) {
      step = -step;
    }

    if (typeofStart == "number") {

      while (step > 0 ? end >= start : end <= start) {
        range.push(start);
        start += step;
      }

    } else if (typeofStart == "string") {

      if (start.length != 1 || end.length != 1) {
        throw TypeError("Only strings with one character are supported.");
      }

      start = start.charCodeAt(0);
      end = end.charCodeAt(0);

      while (step > 0 ? end >= start : end <= start) {
        range.push(String.fromCharCode(start));
        start += step;
      }

    } else {
      throw TypeError("Only string and number types are supported");
    }

    return range;
  }

  /*
   * Return a gaussian function with the given mean and covariance
   */
  fp.normal_pdf = function(mu, sigma) {
    return function(x) {
      return 1.0 / (Math.sqrt(2 * Math.PI) * sigma) * Math.exp(-(x - mu) * (x - mu) / (2 * sigma * sigma));
    };
  };

  fp.plot_function = function(graph_id, chart, f, min, max, name) {
    if (typeof xmin === 'undefined') {
      xmin = 0;
    }
    if (typeof xmax === 'undefined') {
      xmax = 100;
    }
    var xdata = range(xmin, xmax, (xmax - xmin) / 400);

    if (typeof f == 'function') {
      f = [f];
    }

    data = _.zip(f, name).map(function(fn) {
      var func = fn[0], name = fn[1];
      return {
        values: xdata.map(function(x) {return {x: x, y: func(x)};}),
        key: name,
      }
    });

    d3.select(graph_id + ' svg')
      .datum(data)
      .call(chart);
    return chart;
  };

  fp.init_chart = function() {
    var chart  = nv.models.lineChart()
                  .margin({left:100})
                  .useInteractiveGuideline(true)
                  .transitionDuration(350)
                  .showLegend(true)
                  .showYAxis(true)
                  .showXAxis(true);
    chart.xAxis.tickFormat(d3.format('.01f'));
    chart.yAxis.tickFormat(d3.format('.04f'));
    nv.utils.windowResize(function() {chart.update()});
    return chart;
  }

  fp.plot_normal_pdf = function(graph_id, mu, sigma, names, chart) {
    if (typeof chart === 'undefined') {
      chart = fp.init_chart(graph_id);
    }
    if (typeof mu == 'number') {
      mu = [mu];
    }
    if (typeof sigma == 'number') {
      sigma = [sigma];
    }
    if (typeof names == 'string') {
      names = [names];
    }

    if (mu.length != sigma.length) {
      throw TypeError("mu and sigma must have same length.");
    }
    if (typeof names != 'undefined' && mu.length != names.length) {
      throw TypeError("mu and names must have same length.");
    }
    normal_pdfs = _.zip(mu, sigma).map(function(ms) {
      return fp.normal_pdf(ms[0], ms[1]);
    });
    fp.plot_function(graph_id, chart, normal_pdfs, 0, 100, names);
    return chart;
  }

  var matadd = function(A, B) {
    return _.map(_.zip(A, B), function(el) {
      return _.map(_.zip(el[0], el[1]), function(a) {
        return a[0] + a[1];
      });
    });
  }

  var matsub = function(A, B) {
    return _.map(_.zip(A, B), function(el) {
      return _.map(_.zip(el[0], el[1]), function(a) {
        return a[0] - a[1];
      });
    });
  }

  var transpose = function(A) {
    var B = [], i = 0, j = 0;
    for (i = 0; i < A[0].length; i++) {
      B.push([]);
      for (j = 0; j < A.length; j++) {
        B[i].push(A[j][i]);
      }
    }
    return B;
  }

  fp.matmul = function(A, B) {
    var i, j, k;
    var C = [];
    for (i = 0; i < A.length; i++) {
      C.push([])
      for (j = 0; j < B[0].length; j++) {
        C[i].push(0);
        for (k = 0; k < B.length; k++) {
          C[i][j] += A[i][k] * B[k][j];
        }
      }
    }
    return C;
  }

  fp.init_2d = function() {
    var xrange = [0, 100],
      yrange = [-20, 20],
      chart  = function(selection) {

        var data = selection.datum()
        var original_data = data;
        svg = selection.select('svg');
        svg.classed("graphcontainer", true);
        svg.attr('height', 300);
        svg.attr('width', 720);

        var height = parseInt(svg.attr('height'));
        var width = parseInt(svg.attr('width'));
        var rpad = 20;
        var tpad = 20;
        var lpad = 70;
        var bpad = 70;

        var lloc = lpad;
        var rloc = width - rpad;
        var bloc = height - bpad;
        var tloc = tpad;
        var lax_pos = lloc - 60;
        var bax_pos = bloc + 50;
        var midx = (rloc + lloc) / 2;
        var midy = (bloc + tloc) / 2;

        var xscale = d3.scale.linear()
            .domain(xrange)
            .range([lloc, rloc])
            .nice();
        var yscale = d3.scale.linear()
            .domain(yrange)
            .range([bloc, tloc])
            .nice();

        var xAxis = svg.selectAll('.xaxis')
            .data([xscale]);
        var yAxis = svg.selectAll('.yaxis')
            .data([yscale]);

        var enter = xAxis.enter()
            .append("g");
        enter.append("g").classed('axis', true);
        enter.append("g").classed('axislabel', true);

        xAxis.classed('xaxis', true)
            .select(".axis")
              .attr("transform", "translate(0, " + bloc + ")")
              .call(d3.svg.axis().scale(xscale).orient("bottom"))
        xAxis.select('.axislabel')
            .attr("transform", 
                 "translate(" + midx + "," + (bax_pos) + ")")
            .html("<text> Cart Position </text>");

        enter = yAxis.enter()
            .append("g");
        enter.append("g").classed('axis', true);
        enter.append("g").classed('axislabel', true);
        yAxis.classed('yaxis', true)
            .select(".axis")
              .attr("transform", "translate(" + lloc + ",0)")
              .call(d3.svg.axis().scale(yscale).orient("left"))
        yAxis.select('.axislabel')
            .attr("transform", 
                 " rotate(90, " + (lax_pos) + ", " + midy + ") translate(" + (lax_pos) + "," + (midy) + ")")
            .html("<text> Cart Velocity </text>");

        var mapsigma = function(sigma) {
          var xs = xscale(0);
          var ys = yscale(0);
          var ms = [[xscale(xscale(sigma[0][0]) - xs) - xs,
                   yscale(xscale(sigma[0][1]) - xs) - ys],
                  [yscale(xscale(sigma[0][1]) - xs) - ys,
                   yscale(yscale(sigma[1][1]) - ys) - ys]];
          return ms;

        }

        var set_matrix = function(mat) {
          var setfunc = function(selection) {
            var i, j;
            for (i = 0; i < mat.length; i++) {
              for (j = 0; j < mat[i].length; j++) {
                var name = '.mat_' + i + '_' + j;
                selection.select('.mat_' + i + '_' + j)
                  .property('value', mat[i][j]);
              }
            }
          }
          return setfunc;
        }


        /*
         * KF equations
         * d_{t} = d_{t - 1} + v_{t - 1} + 1/2 * u_{t} + 1/2 * eps_v + eps_p
         * v_{t} = v_{t - 1} + u_t + eps_v
         *
         * bar mu_t = A mu_{t-1} + B u_t
         * bar sigma_t = A sigma_{t-1} A^T + R
         */

        var A = [[1, 1], [0, 1]];
        var B = [[0.5], [1]];
        var sigma_v = 1;
        var sigma_p = 0.1;
        var R = [[sigma_p + sigma_v / 4.0, sigma_v / 2], [sigma_v / 2, sigma_v]];
        var C = [[1, 0]];
        var Q = [[3]];

        var display = function() {
          ellipses = svg.selectAll('.sigma')
            .data(data);

          ellipses.enter()
            .append('g').classed('sigma', true)
              .append('ellipse')
                .classed('sigmaellipse', true)
                .attr('style', 'stroke:' + fp.colorwheel[1]);

          ellipses.transition(100)
            .attr("transform",
                  function (x) {
                    var tstr = "translate(" + xscale(x.mu[0]) + "," + yscale(x.mu[1]) + ") ";
                    var ev = fp.eigvecs2x2(mapsigma(x.sigma));
                    var rot = Math.atan2(ev[0][1], ev[0][0]) * 180/Math.PI;
                    var rstr = "rotate(" + rot + ") "
                    return tstr + rstr;
                  });
          ellipses.select('ellipse').transition(100)
            .attr('cx', 0)
            .attr('cy', 0)
            .attr('rx', function(x) {
              return Math.sqrt(fp.eigvals2x2(mapsigma(x.sigma))[0]);
            })
            .attr('ry', function(x) {
              return Math.sqrt(fp.eigvals2x2(mapsigma(x.sigma))[1]);
            })
            .attr('style', function(x, i) {
              return 'stroke:' + fp.colorwheel[(i + 1) % fp.colorwheel.length];
            })



          selection.select('.mean')
            .call(set_matrix(data[0].mu));
          selection.select('.covariance')
            .call(set_matrix(data[0].sigma));

          selection.select('.ivec')
            .call(set_matrix(fp.matmul(fp.matinv2x2(data[0].sigma), data[0].mu)));
          selection.select('.imat')
            .call(set_matrix(fp.matinv2x2(data[0].sigma)));

          selection.select('.movementnoise')
            .call(set_matrix(R));
          selection.select('.measurementnoise')
            .call(set_matrix(Q));
        }

        display();


        var apply_force = function(x) {
          R = read_matrix(selection.select('.movementnoise'));
          Q = read_matrix(selection.select('.measurementnoise'));

          data = data.map(function(y) {
            return {
              key: y.key,
              mu: matadd(fp.matmul(A, y.mu), fp.matmul(B, [[x]])),
              sigma: matadd(fp.matmul(A, fp.matmul(y.sigma, transpose(A))), R),
            }
          });
          display();
        }
        var shift_left = function() {
          apply_force(-5);
        }
        var shift_right = function() {
          apply_force(5);
        }
        var propagate = function() {
          apply_force(0);
        }
        var reset = function() {
          data = original_data;
          display();
        }

        var mat_entry = function(i, j) {
          return '.mat_' + i + '_' + j;
        }

        var read_matrix = function(selection) {
          var i = 0, j = 0;
          var mat = [];
          while (!selection.select(mat_entry(i, 0)).empty()) {
            mat.push([]);
            j = 0;
            while (!selection.select(mat_entry(i, j)).empty()) {
              mat[i].push(parseFloat(selection.select(mat_entry(i, j)).property('value')));
              j ++;
            }
            i ++;
          }
          return mat;
        }

        var read_matrices = function() {
        }

        // selection.selectAll('input').on('change', read_matrices);

        selection.select('.rightshift').on('click', shift_right);
        selection.select('.leftshift').on('click', shift_left);
        selection.select('.continue').on('click', propagate);
        selection.select('.reset').on('click', reset);

        svg.on('mousemove', function() {
          var scoord = d3.mouse(this);
          var x = xscale.invert(scoord[0])
          var indicator = svg.selectAll('.indicator')
            .data([x]);

          indicator.enter()
            .append('line')
            .classed('indicator', true)
            .attr('y1', bloc)
            .attr('y2', tloc);

          indicator.attr('x1', xscale(x));
          indicator.attr('x2', xscale(x));
        });

        svg.on('click', function () {
          var scoord = d3.mouse(this);
          var x = xscale.invert(scoord[0])
          var y = yscale.invert(scoord[0])
          var measurement = [[xscale.invert(scoord[0])]];
          R = read_matrix(selection.select('.movementnoise'));
          Q = read_matrix(selection.select('.measurementnoise'));
          data = data.map(function(x) {
            // K = sigma * C^T * (C sigma C^T + Q)^-1
            var K = fp.matmul(fp.matmul(x.sigma, transpose(C)),
                              fp.matinv2x2(matadd(fp.matmul(fp.matmul(C, x.sigma), transpose(C)), Q)));
            // mu = mu_{t-1} + K (z - C mu_{t-1}
            var mu = matadd(x.mu, fp.matmul(K, matsub(measurement, fp.matmul(C, x.mu))));
            var sigma = fp.matmul(matsub([[1, 0], [0, 1]], fp.matmul(K, C)), x.sigma);
            return {
              key: x.key,
              mu: mu,
              sigma: sigma,
            }
          });
          display();
        });
      }

    return chart;
  }

  fp.eigvals2x2 = function(mat) {
    var a = mat[0][0],
        b = mat[0][1],
        c = mat[1][0],
        d = mat[1][1],
        offset, rec, eigenvalues;
        diff = (a + d) * (a + d) - 4 * (a * d - b * c);
        if (diff > 0)
        {
          offset = (a + d) / 2;
          rec = Math.sqrt(diff) / 2;
        }
        else
        {
          offset = (a + d) / 2;
          rec = 0;
        }
        eigenvalues = [
            offset + rec,
            offset - rec
        ];
    return eigenvalues;
  }

  fp.matinv2x2 = function(mat) {
    if (mat.length == 1) {
        var inv = [[1.0 / mat[0][0]]];
    } else {
        var det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0],
            inv = [
                [mat[1][1] / det, - mat[1][0] / det],
                [- mat[0][1] / det, mat[0][0] / det]
            ];
    }
    return inv;
  }

  fp.eigvecs2x2 = function(mat) {
    /* 
     * There is something really wrong with this eigenvector
     * calculation, I need to fix it.
     **/
    var ev = fp.eigvals2x2(mat),
      eigenvecs, i;

    if (Math.abs(mat[0][1]) > 0) {
      eigenvecs = [
            [mat[0][1], ev[0] - mat[0][0]],
            [mat[0][1], ev[1] - mat[0][0]]
        ];
    } else if (Math.abs(mat[1][0]) > 0) {
      eigenvecs = [
          [ev[0] - mat[1][1], ev[1][0]],
          [ev[1] - mat[1][1], ev[1][0]]
        ];
    } else {
      // Diagonal, so eigenvectors are identity, or rotated identity
      if (Math.abs(mat[0][0] - ev[0]) < Math.abs(mat[0][0] - ev[1])) {
        eigenvecs = [
          [1, 0],
          [0, 1]
        ]
      } else {
        eigenvecs = [
            [0, 1],
            [1, 0]
        ];
      }
    }

    for (i = 0; i < 2; i++) {
        var tmp = Math.sqrt(Math.pow(eigenvecs[i][0], 2)
                            + Math.pow(eigenvecs[i][1], 2))
        eigenvecs[i][0] /= tmp;
        eigenvecs[i][1] /= tmp;
    }

    return eigenvecs;
  }

  fp.plot_2d = function(selector, chart, mean, covariance, names) {

    var pos_scale = d3.scale.linear()
              .domain(pos_range)
              .range([])
              .nice(),
        vel_scale = d3.scale.linear()
              .domain(pos_range)
              .range([])
              .nice(),
        xAxis = d3.svg.axis()
          .scale(pos_scale)
          .orient("bottom"),
        yAxis = d3.svg.axis()
          .scale(vel_scale)
          .orient("left");
        mmean = [xscale(mean[0])]
    if (typeof chart === 'undefined') {
      chart = fp.init_2d(selector);
    }
  }

  /* Nonlinear comparison stuff */
  fp.init_nonlinear_comparison = function() {
    var xrange = [-20, 20],
      yrange = [-20, 20],
      chart = function(selection) {
        var data = selection.datum(),
          original_data = data,
          svg = selection.select('svg');
        selection.classed("graphcontainer", true);
        selection.classed("square", true);
        svg.classed("graphcontainer", true);
        svg.classed("square", true);
        svg.attr('height', 500);
        svg.attr('width', 500);

        var height = parseInt(svg.attr('height')),
          width = parseInt(svg.attr('width')),
          rpad = 20,
          tpad = 20,
          lpad = 70,
          bpad = 70,
          lloc = lpad,
          rloc = width - rpad,
          bloc = height - bpad,
          tloc = tpad,
          lax_pos = lloc - 60,
          bax_pos = bloc + 50,
          midx = (rloc + lloc) / 2,
          midy = (bloc + tloc) / 2,
          Q = [[1, 0], [0, 1]],
          xscale = d3.scale.linear()
            .domain(xrange)
            .range([lloc, rloc])
            .nice(),
          yscale = d3.scale.linear()
            .domain(yrange)
            .range([bloc, tloc])
            .nice(),
          xAxis = svg.selectAll('.xaxis')
            .data([xscale]),
          yAxis = svg.selectAll('.yaxis')
            .data([yscale]),
          enter = xAxis.enter()
            .append("g");
        enter.append("g").classed('axis', true);
        enter.append("g").classed('axislabel', true);

        xAxis.classed('xaxis', true)
            .select(".axis")
              .attr("transform", "translate(0, " + bloc + ")")
              .call(d3.svg.axis().scale(xscale).orient("bottom"))
        xAxis.select('.axislabel')
            .attr("transform", 
                 "translate(" + midx + "," + (bax_pos) + ")")
            .html("<text> Robot X Position </text>");

        enter = yAxis.enter()
            .append("g");
        enter.append("g").classed('axis', true);
        enter.append("g").classed('axislabel', true);
        yAxis.classed('yaxis', true)
            .select(".axis")
              .attr("transform", "translate(" + lloc + ",0)")
              .call(d3.svg.axis().scale(yscale).orient("left"))
        yAxis.select('.axislabel')
            .attr("transform", 
                 " rotate(90, " + (lax_pos) + ", " + midy + ") translate(" + (lax_pos) + "," + (midy) + ")")
            .html("<text> Robot Y Position </text>");
        /*
         * Try drawing lines from mean to boundary of ellipse at
         * the mean
         *
         * direction of line is phy = theta - rot
         * length of line is the distance to the boundary of the ellipse
         * in that direction....
         * x^2 / rx^2 + y^2 / ry^2 = 1
         * cos(theta)^2 / rx^2 + 
         */

        var c_xscale_inv = function(x) {
          var xs = xscale.invert(0);
          var ret =  xscale.invert(x) - xs;
          return ret;
        }
        var c_yscale_inv = function(y) {
          var ys = yscale.invert(0);
          return -(yscale.invert(y) - ys);
        }
        var c_xscale = function(x) {
          var xs = xscale(0);
          return xscale(x) - xs
        };
        var c_yscale = function(y) {
          var ys = yscale(0);
          return -(yscale(y) - ys);
        };
        var extract_mapsigma = function(sigma) {
          var xs = xscale(0),
            ys = yscale(0),
            ms = [[c_xscale(c_xscale(sigma[0][0])),
                   c_yscale(c_xscale(sigma[0][1]))],
                  [c_xscale(c_yscale(sigma[1][0])),
                   c_yscale(c_yscale(sigma[1][1]))]];
          return ms;

        }

        var display = function() {
          robots = svg.selectAll('.sigma')
            .data(data);

          var new_g = robots.enter()
            .append('g').classed('sigma', true);

          new_g.append('line')
                .classed('plussigma', true)
                .classed('sigmaline', true)
                .attr('style', 'stroke:' + fp.colorwheel[0]);
          new_g.append('line')
                .classed('minussigma', true)
                .classed('sigmaline', true)
                .attr('style', 'stroke:' + fp.colorwheel[0]);
          new_g.append('line')
                .classed('mean_theta', true)
                .classed('sigmaline', true)
                .attr('style', 'stroke:' + fp.colorwheel[1]);
          new_g.append('ellipse')
                .classed('sigmaellipse', true)
                .attr('style', 'stroke:' + fp.colorwheel[1]);

          robots.transition(100)
            .attr("transform",
                  function (x) {
                    var tstr = "translate(" + xscale(x.mu[0]) + "," + yscale(x.mu[1]) + ") ";
                    var ev = fp.eigvecs2x2(extract_mapsigma(x.sigma));
                    var rot = Math.atan2(ev[0][1], ev[0][0]) * 180/Math.PI;
                    var rstr = "rotate(" + rot + ") "
                    x.rotation = rot * Math.PI / 180;
                    return tstr + rstr;
                  });

          robots.select('ellipse').transition(100)
            .attr('cx', 0)
            .attr('cy', 0)
            .attr('rx', function(x) {
              return Math.sqrt(fp.eigvals2x2(extract_mapsigma(x.sigma))[0]);
            })
            .attr('ry', function(x) {
              return Math.sqrt(fp.eigvals2x2(extract_mapsigma(x.sigma))[1]);
            })
            .attr('style', function(x, i) {
              return 'stroke:' + fp.colorwheel[(i + 1) % fp.colorwheel.length];
            });

          var get_phi = function(x, start_phi) {
            /*
             * Map to raw (true) rotation
             */
            var mapped_rot = x.rotation,
              raw_cos = c_xscale_inv(Math.cos(mapped_rot)),
              raw_sin = c_yscale_inv(Math.sin(mapped_rot)),
              raw_rot = Math.atan2(raw_sin, raw_cos),
            /*
             * Add the desired angle and map back to displayed rotation
             */
              x = c_xscale(Math.cos(-start_phi - raw_rot)),
              y = c_yscale(Math.sin(-start_phi - raw_rot)),
              phi = Math.atan2(y, x);
              return phi;
          }
          var get_radii = function(x) {
            var ev1 = fp.eigvals2x2(extract_mapsigma(x.sigma)),
              rx = Math.sqrt(ev1[0]),
              ry = Math.sqrt(ev1[1]);
            return [rx, ry];
          }
          robots.select('.minussigma').transition(100)
            .attr('x1', 0)
            .attr('y1', 0)
            .attr('x2', function(x) {
              var phi = get_phi(x, x.mu[2][0] - Math.sqrt(x.sigma[2][2])),
                rx = get_radii(x)[0],
                phix = rx * Math.cos(phi);
              return phix;
            })
            .attr('y2', function(x) {
              var phi = get_phi(x, x.mu[2][0] - Math.sqrt(x.sigma[2][2])),
                ry = get_radii(x)[1],
                phiy = ry * Math.sin(phi);
              return phiy;
            });
          robots.select('.plussigma').transition(100)
            .attr('x1', 0)
            .attr('y1', 0)
            .attr('x2', function(x) {
              var phi = get_phi(x, x.mu[2][0] + Math.sqrt(x.sigma[2][2])),
                rx = get_radii(x)[0],
                phix = rx * Math.cos(phi);
              return phix;
            })
            .attr('y2', function(x) {
              var phi = get_phi(x, x.mu[2][0] + Math.sqrt(x.sigma[2][2])),
                ry = get_radii(x)[1],
                phiy = ry * Math.sin(phi);
              return phiy;
            });
          robots.select('.mean_theta').transition(100)
            .attr('x1', 0)
            .attr('y1', 0)
            .attr('x2', function(x) {
              var phi = get_phi(x, x.mu[2][0]),
                rx = get_radii(x)[0],
                phix = rx * Math.cos(phi);
              return phix;
            })
            .attr('y2', function(x) {
              var phi = get_phi(x, x.mu[2][0]);
                ry = get_radii(x)[1],
                phiy = ry * Math.sin(phi);
              return phiy;
            });
        };

        display();

        svg.on('mousemove', function() {
          var scoord = d3.mouse(this);
          var x = xscale.invert(scoord[0]),
            y = yscale.invert(scoord[1]),
            xindicator = svg.selectAll('.xindicator')
              .data([x]),
            yindicator = svg.selectAll('.yindicator')
              .data([y]);
          x = Math.min(xrange[1], Math.max(xrange[0], x));
          y = Math.min(yrange[1], Math.max(yrange[0], y));

          xindicator.enter()
            .append('line')
            .classed('indicator', true)
            .classed('xindicator', true)
            .attr('y1', bloc)
            .attr('y2', tloc);

          yindicator.enter()
            .append('line')
            .classed('indicator', true)
            .classed('yindicator', true)
            .attr('x1', lloc)
            .attr('x2', rloc);

          xindicator.attr('x1', xscale(x));
          xindicator.attr('x2', xscale(x));
          yindicator.attr('y1', yscale(y));
          yindicator.attr('y2', yscale(y));

          ellipses = svg.selectAll('.Qellipse')
            .data([{'x': x, 
                    'y': y, 
                    'sigma': Q
            }]);
          
          ellipses.enter()
            .append('g').classed('Qellipse', true)
              .append('ellipse')
                .classed('sigmaellipse', true)
                .classed('Qellipse', true);

          ellipses.attr("transform", function (x) {
              var tstr = "translate(" + xscale(x.x) + "," + yscale(x.y) + ") ",
                ev = fp.eigvecs2x2(extract_mapsigma(x.sigma)),
                rot = Math.atan2(ev[0][1], ev[0][0]) * 180/Math.PI,
                rstr = "rotate(" + rot + ") ";
              return tstr + rstr;
            });

          ellipses.select('ellipse')
            .attr('cx', 0)
            .attr('cy', 0)
            .attr('rx', function(x) {
              return Math.sqrt(fp.eigvals2x2(extract_mapsigma(Q))[0]);
            })
            .attr('ry', function(x) {
              return Math.sqrt(fp.eigvals2x2(extract_mapsigma(Q))[1]);
            });
        });

        svg.on('click', function() {
          var scoord = d3.mouse(this);
          var x = xscale.invert(scoord[0]),
            y = yscale.invert(scoord[1]),
            xindicator = svg.selectAll('.xindicator')
              .data([x]),
            yindicator = svg.selectAll('.yindicator')
              .data([y]);
          x = Math.min(xrange[1], Math.max(xrange[0], x));
          y = Math.min(yrange[1], Math.max(yrange[0], y));

          var alldata = {
            robots: data,
            Q: Q,
            x: x,
            y: y
          };
          
          d3.json("/measurementupdate/ukf")
            .header("Content-Type", "application/json")
            .post(JSON.stringify(alldata), function(error, curdata) {
              if (error) {
                return;
              }
              data = curdata.robots;
              console.log(data)
              display();
            });
        });

        var movement_update = function(leftm, rightm) {
          var alldata = {
            robots: data,
            leftwheel: leftm,
            rightwheel: rightm
          }

          d3.json("/movementupdate/ukf")
            .header("Content-Type", "application/json")
            .post(JSON.stringify(alldata), function(error, curdata) {
              if (error) {
                return;
              }
              data = curdata.robots;
              display();
            });
        }

        var movebutton = selection.select('.movebutton'),
          leftmove = selection.select('.leftwheel'),
          rightmove = selection.select('.rightwheel');
        console.log(movebutton);
          // Update button events
        selection.select('.movebutton').on('click', function() {
            var leftm = leftmove.property('value'),
              rightm = rightmove.property('value');
            console.log("Clicked");
            movement_update(leftm, rightm);
            display();
          });


      };
      return chart;
  }

  return fp;
}();

gaussianinterface = function(selector, options) {
  var default_options = {
    movement_buttons: false,
    noise_selection: false,
    mu_setting: false,
    sigma_setting: false,
    invec_setting: false,
    inmat_setting: false,
    sensor_buttons: false,
    default_values: {
      mu: [40, 40, 40],
      sigma: [1, 4, 20],
      names: ['localized', 'moderate', 'delocalized'],
    },
  }
  var opts = $.extend(default_options, options);
  var template = Handlebars.compile($('#graph_template').html());
  var contents = template(opts);
  var gi = {mu: opts.default_values.mu, sigma: opts.default_values.sigma,
    names: opts.default_values.names};
  $(selector).html(contents);

  gi.redisplay = function() {
    filter_plotting.plot_normal_pdf(selector, gi.mu, gi.sigma, gi.names, gi.chart);
    gi.show_values();
  };

  gi.show_values = function() {
    var mu_d = gi.el.select(".mean");
    var sigma_d = gi.el.select('.covariance');
    var xi_d = gi.el.select('.information_vec');
    var omega_d = gi.el.select('.information_mat');
    if (!mu_d.empty()) {
      mu_d.property('value', gi.mu[0]);
    }
    if (!sigma_d.empty()) {
      sigma_d.property('value', gi.sigma[0]);
    }
    if (!xi_d.empty()) {
      xi_d.property('value', gi.mu / gi.sigma[0]);
    }
    if (!omega_d.empty()) {
      omega_d.property('value', 1 / gi.sigma[0]);
    }
  }

  gi.getnoise = function() {
    var tmp = d3.select(selector + ' .noise');
    if (!tmp.empty()) {
      var noise = tmp.property('value');
      if (noise < 0) {
        noise = -noise;
      }
      tmp.property('value', noise);
      return noise;
    }
    return 0;
  }

  gi.shift = function(dmu) {
    if (typeof dmu === 'undefined') {
      dmu = -5;
    }
    var noise = gi.getnoise();
    gi.sigma = gi.sigma.map(function(sigma) {return Math.sqrt(sigma * sigma + noise * noise);});
    gi.mu = gi.mu.map(function(mu) {return mu + dmu;});
    gi.redisplay();
  }

  gi.shift_right = function(dmu) {
    if (typeof dmu === 'undefined') {
      dmu = 5;
    }
    gi.shift(dmu);
  }

  gi.shift_left = function(dmu) {
    if (typeof dmu === 'undefined') {
      dmu = 5;
    }
    gi.shift(-dmu);
  }

  gi.reset = function() {
    gi = $.extend(gi, opts.default_values);
    gi.redisplay();
  }

  
  gi.chart = filter_plotting.plot_normal_pdf(selector, gi.mu, gi.sigma, gi.names);
  gi.chart.xAxis.axisLabel('Distance along track: d (meters)');
  gi.chart.yAxis.axisLabel('Probability density: p(d)')
  gi.selector = selector;
  gi.el = d3.select(selector);
  gi.show_values();

  gi.el.select('.leftshift')
    .on('click', gi.shift_left);
  gi.el.select('.rightshift')
    .on('click', gi.shift_right);
  gi.el.select('.reset')
    .on('click', gi.reset);
}
