
FluxPalette <- function(n,type=1:3){
	colorRampPalette(
		switch(type[1],
			c("#053061",
			"#2166ac",
			"#4393c3",
			"#92c5de",
			"#d1e5f0",
			"#f7f7f7",
			"#fddbc7",
			"#f4a582",
			"#d6604d",
			"#b2182b",
			"#67001f")
			,c("#313695",
			"#4575B4",
			"#74ADD1",
			"#ABD9E9",
			"#E0F3F8",
			"#FFFFBF",
			"#FEE090",
			"#FDAE61",
			"#F46D43",
			"#D73027",
			"#A50026")
			,c("#ffffd9",
			"#edf8b1",
			"#c7e9b4",
			"#7fcdbb",
			"#41b6c4",
			"#1d91c0",
			"#225ea8",
			"#253494",
			"#081d58")
			)
	,bias=1,interpolate="linear",space="Lab")(n)
}

ConcPalette <- colorRampPalette(c(
	"#006837",
	"#1a9850",
	"#66bd63",
	"#a6d96a",
	"#d9ef8b",
	"#ffffbf",
	"#fee08b",
	"#fdae61",
	"#f46d43",
	"#d73027",
	"#a50026"),bias=1,interpolate="linear",space="Lab")

