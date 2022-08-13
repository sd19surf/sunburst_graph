

    var hourScale = d3.scaleLinear()
        .range([0,345])
        .domain([0,23]);
    var radians = 0.0174532925;
    var currentDate = new Date();
    var i_year = currentDate.getFullYear();
    var i_month = currentDate.getMonth() + 1;
    var i_day = currentDate.getDate();
    var i_lat = 33.973057;
    var i_long = -80.472778;


    var solarData = json_main(i_year,i_month,i_day,i_long,i_lat,0).events;
    console.log(solarData);
    var d3Data = constructSolarData(solarData);
    console.log(d3Data);
    var svg = d3.select("svg")
        .append("g")
        .attr("transform", "translate(150,150)");
    
    // clock tickers
    var face = svg.append('g')
        .attr('id','clock-face');

    face.selectAll('.hour-tick')
        .data(d3.range(0,24)).enter()
        .append('line')
        .attr('class', 'hour-tick')
        .attr('x1',0)
        .attr('x2',0)
        .attr('y1',120)
        .attr('y2',128)
        .attr('transform',function(d){
            return 'rotate(' + hourScale(d) + ')';
        });

    face.selectAll('.hour-label')
        .data(d3.range(0,24,1))
        .enter()
        .append('text')
        .attr('class', 'hour-label')
        .attr('text-anchor','middle')
        .attr('x',function(d){
            return 135*Math.sin(hourScale(d)*radians);
        })
        .attr('y',function(d){
            return -135*Math.cos(hourScale(d)*radians)+4;
        })
        .text(function(d){
            return d+"Z";
        });

    // The gray background for the chart
    var backgroundPath = d3.arc()
        .innerRadius(110)
        .outerRadius(120)
        .startAngle(0)
        .endAngle((Math.PI/180) * 360);

    var background = svg.append("path")
        .attr("class","backdrop")
        .attr("d", backgroundPath)
        .attr("fill", "#000")

    // Appending the sections paths as arcs
// Container for the arcs
    var sectors = svg.append("g");

    sectors.selectAll("path")
        .data(d3Data)
        .attr("class","arc")
        .enter().append("path")
        .each(arcFunction);




// Function called for each path appended to increase scale and iterate.
function arcFunction(d, i) {
var arc = d3.arc()
        .innerRadius(110)
        .outerRadius(120)
        .startAngle(d.angles.startAngle)
        .endAngle(d.angles.endAngle)
 


return d3.select(this)
    .style("fill", d.color)
    .attr("d", arc);
} 

function constructSolarData(rawData){
    newData = [];
    newData.push({"name":"MR","color": "silver","angles":timeToRadians(solarData[0].time,solarData[1].time)});
     newData.push({"name":"BMAT", "color":"grey", "angles":timeToRadians(solarData[8].time,solarData[6].time)});
    newData.push({"name":"BMNT", "color":"lightgrey", "angles":timeToRadians(solarData[6].time,solarData[4].time)});
    newData.push({"name":"BMCT", "color":"lightyellow", "angles":timeToRadians(solarData[4].time,solarData[2].time)});
    newData.push({"name":"EENT", "color":"lightyellow", "angles":timeToRadians(solarData[3].time,solarData[5].time)});
    newData.push({"name":"EECT", "color":"lightgrey", "angles":timeToRadians(solarData[5].time,solarData[7].time)});
    newData.push({"name":"EEAT", "color":"grey", "angles":timeToRadians(solarData[7].time,solarData[9].time)});
    newData.push({"name":"SR", "color":"yellow", "angles":timeToRadians(solarData[2].time,solarData[3].time)});
    return newData;

}
function timeToRadians(startTime="",endTime=""){

    dateStart = new Date(startTime);
    dateEnd = new Date(endTime);
    
    dateEndDecimal = (dateEnd.getUTCHours() + dateEnd.getUTCMinutes() / 60);
    dateStartDecimal = (dateStart.getUTCHours() + dateStart.getUTCMinutes() / 60);
    radianStartHours = ((360 / 24) * dateStartDecimal)*(Math.PI/180);
    radianEndHours = ((360 / 24) *dateEndDecimal)*(Math.PI/180);
    if (radianStartHours < radianEndHours){
    return {'startAngle':radianStartHours,'endAngle':radianEndHours};
    }else{
        radianEndHours = ((360/24)* 23.9)*(Math.PI/180);
    return {'startAngle':radianStartHours,'endAngle':radianEndHours};
    }
}
