// d3 = require('d3.min.js');

data = d3.json('tutorial_example.json');

var margin = {top: 30, right: 20, bottom: 30, left: 50},
width = 600 - margin.left - margin.right,
height = 270 - margin.top - margin.bottom;


// Set the ranges
var x = d3.scale.linear().range([0, width]);
var y = d3.scale.linear().range([height, 0]);

// Define the axes
var xAxis = d3.svg.axis().scale(x)
.orient("bottom").ticks(5);

var yAxis = d3.svg.axis().scale(y)
.orient("left").ticks(5);

// Define the line
var valueline = d3.svg.line()
.x(function(d) { return x(d.date); })
.y(function(d) { return y(d.close); });

// Adds the svg canvas
var svg = d3.select("body")
.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
.append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");

// Get the data
data = d3.json('tutorial_example.json', function(error, data){
    console.log(data);
    for(i = 0; i < data.electrolytes[0][1].concentrations.length; i++){
        data.date = data.electrolytes[0][1].nodes;
        data.close = data.electrolytes[0][1].concentrations[i];
        data2=[]
        for(i = 0; i < data.date.length; i++){
            data2[i] = {date: data.date[i], close: data.close[i]}
        }


        // Scale the range of the data
        x.domain([0, d3.max(data2, function(d) { return d.date; })]);
        y.domain([0, d3.max(data2, function(d) { return d.close; })]);


        // Add the valueline path.
        svg.append("path")
            .attr("class", "line")
            .attr("d", valueline(data2))

        // Add the X Axis
        svg.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + height + ")")
            .call(xAxis);

        // Add the Y Axis
        svg.append("g")
            .attr("class", "y axis")
            .call(yAxis);

        }});
