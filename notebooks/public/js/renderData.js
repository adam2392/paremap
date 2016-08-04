// data files
var faces_csv = './data/faces.csv';
var roi_csv = './data/roi.csv';
var vertices_csv = './data/vertices.csv';

// initialize size params
var m_width = $("#scatterplots").width();
var m_height = $("#scatterplots").height();
var width = 938;
// var mapRatio = 0.6;
var height = 500;

// var margin = {top:30, right:20, bottom:30, left:50},
// 	width = 600 - margin.left - margin.right,
// 	height = 270 - margin.top - margin.bottom;


// create the svg_pol root inside our div html
var svg = d3.select("#scatterplots").append("svg")
	.attr("width", m_width)
	.attr("height", m_width*height/width)
	.attr("preserveAspectRatio", "xMidYMid")
	.attr("viewBox", "0 0 " + width + " " + height)

// Optional: Add background behind the displayed map
svg.append("rect")
    .attr("class", "background")
    .attr("fill", "white")
    .attr("width", width)
    .attr("height", height);

var g = svg.append("g");


d3.queue()
	.defer(d3.csv, faces_csv)
	.defer(d3.csv, roi_csv)
	.defer(d3.csv, vertices_csv)
	.await(ready);

/* MAIN FUNCTION AREA */
function ready(error, faces, roi, vertices) {
	if (error) console.log(error);

	// get x,y,z coordinates of each vertex
	var x_point = vertices.map(function(value,index) { 
		return value['x']; });
	var y_point = vertices.map(function(value,index) { 
		return value['y']; });
	var z_point = vertices.map(function(value,index) { 
		return value['z']; });

	// console.log(x_point)
	// console.log(y_point)
	// console.log(z_point)
	// console.log(vertices)

	g.selectAll("circle")
		.data(vertices).enter()
		.append("circle")
		.attr("cx", function(d) { 
				// console.log(d)
				return d['x'];
			})
		.attr("cy", function(d) {
				return d['y'];
			})
		.attr("cz", function(d) {
				return d['z'];
			})
		.attr("r", "2px")
		.attr("fill", "gray")
		.attr("transform", "translate(" + 200 + "," + 100 + ")")
}
