var faces_csv = "./data/faces.csv"
var roi_csv = "./data/roi.csv"
var vertices_csv = "./data/vertices.csv"

// 0. Define parameters
var m_width = $("#scatterplots").width()
var width = 938;
var height = 500;

// 1. Create SVG
var svg = d3.select("scatterplots").append("svg")
	.attr("width", width)
	.attr("height", height)
	.attr("preserveAspectRatio", "xMidYMid")
	.attr("viewBox", "0 0 " + width + " " + height)

var g = svg.append("g")

// 2. Queue up the data
d3.queue()
	.defer(d3.csv, faces_csv)
	.defer(d3.csv, roi_csv)
	.defer(d3.csv, vertices_csv)
	.await(ready)

// Main Function
function ready(error, faces, roi, vertices) {
	if (error) throw error;

	
}