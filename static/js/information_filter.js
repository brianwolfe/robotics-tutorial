
gaussianinterface('#information_interface', 
                  {
                    movement_buttons: true,
                    noise_selection:true,
                    mu_setting: true,
                    sigma_setting: true,
                    invec_setting: true,
                    inmat_setting: true,
                    default_values: {
                      mu: [40],
                      sigma: [4],
                      names: ['probability'],
                    }
                  });


var tmp = 1;
var chart = filter_plotting.init_2d()
var template = Handlebars.compile($('#twod_template').html());

Handlebars.registerPartial('vector', $('#vector-partial').html());
Handlebars.registerPartial('matrix', $('#matrix-partial').html());
Handlebars.registerPartial('mat11', $('#mat11-partial').html());

var contents = template({});
d3.select("#twod").html(contents);

d3.select("#twod")
  .datum([
    {
      key: "initial",
      mu: [[40], [0]],
      sigma: [[10, 0], [0, 20]],
    }
  ])
  .call(chart);


