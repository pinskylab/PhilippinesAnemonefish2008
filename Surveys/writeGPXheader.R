# Write the header to a GPX file

writeGPXheader = function(filename){
	date = date()

	header=c(
	"<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"yes\"?>",
	"<gpx",
 	" version=\"1.0\"",
	" creator=\"R 2.7.1 writeGPXheader.R\"", 
	" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"",
	" xmlns=\"http://www.topografix.com/GPX/1/0\"",
	" xmlns:topografix=\"http://www.topografix.com/GPX/Private/TopoGrafix/0/2\"",
	" xsi:schemaLocation=\"http://www.topografix.com/GPX/1/0 http://www.topografix.com/GPX/1/0/gpx.xsd http://www.topografix.com/GPX/Private/TopoGrafix/0/2 http://www.topografix.com/GPX/Private/TopoGrafix/0/2/topografix.xsd\">",
	" <author>Malin Pinsky</author>",
	" <email>mpinsky@stanford.edu</email>",
	" <url>an_url</url>",
	" <urlname>a_url_name</urlname>",
	paste("<time>",date,"</time>", sep="")
	)

	writeLines(header, filename, sep="\n")

}