// set plot size
var margin = {top: 30, right: 20, bottom: 30, left: 100},
width = 600 - margin.left - margin.right,
height = 270 - margin.top - margin.bottom;


// Set the ranges
var x = d3.scale.linear().range([0, width]);
var y = d3.scale.linear().range([height, 0]);

// Define the axes
var xAxis = d3.svg.axis()
                  .scale(x)
                  .orient("bottom")
                  .ticks(5);

var yAxis = d3.svg.axis()
                  .scale(y)
                  .orient("left")
                  .ticks(5);

// Define the line
var valueline = d3.svg.line()
                      .x(function(d) { return x(d[0]); })
                      .y(function(d) { return y(d[1]); });

// Adds the svg canvas
var svg = d3.select("body")
            .append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// Set up color options
var p=d3.scale.category10();

// Function for setting axis attributes
var make_axes = function(data2){
    // Scale the range of the data
    x.domain([0, d3.max(data2, function(d) { return d[0]; })]);
    y.domain([0, 2*d3.max(data2, function(d) { return d[1]; })]);

    // Add the X Axis
    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

    // Add the Y Axis
    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);

    //add xlabel
    svg.append("text")
       .attr("class", "x label")
       .attr("text-anchor", "end")
       .attr("x", width)
       .attr("y", height - 6)
       .text("length (m)");

    //add ylabel
    svg.append("text")
       .attr("class", "y label")
       .attr("text-anchor", "end")
       .attr("y", -50)
       .attr("dy", ".75em")
       .attr("transform", "rotate(-90)")
       .text("concentration (M)");
}

// Plot button info
var frame = 0
var lines = []
var saved_data = null
var plot = function(){
    frame = 0
    var file = document.getElementById("file_input")
                       .value
                       .replace("C:\\fakepath\\", "");

    plotter = d3.json(file, function(error, data){
        saved_data = data;
        for(i = 0; i < data.electrolytes[frame][1].concentrations.length; i++){
            if(i === 0){make_axes(d3.zip(
                data.electrolytes[frame][1].nodes,
                data.electrolytes[frame][1].concentrations[i]
            ))};

            // Add the valueline path.
            lines[i] = svg.append("path")
                            .attr("class", "line")
                            .attr("d", valueline(d3.zip(
                                data.electrolytes[frame][1].nodes,
                                data.electrolytes[frame][1].concentrations[i]
                            )))
                            .attr("stroke", p(i%10))

        }});
    };

var next = function(time){
    frame++;
    update_lines(frame, time);
};

var prev = function(time){
    frame--;
    update_lines(frame, time);
}

var update_lines = function(frame, time){
    for(i = 0; i < saved_data.electrolytes[frame][1].concentrations.length; i++){
        lines[i].transition().duration(time)
                .attr("d", valueline(d3.zip(
                    saved_data.electrolytes[frame][1].nodes,
                    saved_data.electrolytes[frame][1].concentrations[i])))}};

var play = function(){
    while(frame<saved_data.electrolytes.length){
        next(1000)
    }
}
