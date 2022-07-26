function categories(container) {
  var categories_json;
  container.add_button = () => {
    node = { name: "some_element" };
    categories_json.children.push(node);
    set_hierarchy(json2hierarchy(categories_json));
    reset_ui();
  };
  let width = window.innerWidth,
    height = window.innerHeight,
    radius = Math.min(width, height) / 2 - 10,
    x = d3.scaleLinear().range([0, 2 * Math.PI]),
    y = d3.scaleSqrt().range([0, radius]),
    color = d3.scaleOrdinal(d3.schemeAccent),
    arc = d3
      .arc()
      .startAngle(d => {
        return (d.startAngle = Math.max(0, Math.min(2 * Math.PI, x(d.x0))));
      })
      .endAngle(d => {
        d.addStrokeSpace = false;
        let children_shown = 0;
        if (d.parent && d.parent.children)
          d.parent.children.map(child => {
            if (child.shown) children_shown++;
          });
        if (d != rootNode && children_shown > 1) d.addStrokeSpace = true;
        return (d.endAngle = Math.max(
          0,
          Math.min(2 * Math.PI, x(d.x1) - (d.addStrokeSpace ? 0.007 : 0))
        ));
      })
      .innerRadius(d => {
        return (d.innerR = Math.max(0, y(d.y0) + (d.y0 ? 2 : 0)));
      }) // .svgaon.stroke-width()
      .outerRadius(d => {
        return (d.outerR = Math.max(0, y(d.y1)));
      }),
    currScale = 1.0,
    currX = width / 2,
    currY = height / 2,
    svg = d3
      .select(container)
      .append("svg")
      .attr("width", width)
      .attr("height", height)
      .append("g")
      .attr("transform", "translate(" + currX + "," + currY + ")") // center svg
      .call(
        d3
          .drag()
          .on("drag", dragged)
          .on("start", dragstart)
          .on("end", dragend)
      ),
    formatNumber = d3.format(",d"),
    tooltipShowTime = 200,
    tooltipHideTime = 500,
    tooltip = d3
      .select(container)
      .append("div")
      .attr("class", "tooltip")
      .style("opacity", 0),
    show_tooltips = true,
    tweenTransitionTime = 300,
    last_node_clicked = null;
  click = d => {
    if (d === last_node_clicked) return;
    if (last_node_clicked || d.parent) click_tween(d);
    last_node_clicked = d;
  };
  zoom = d3
    .zoom()
    .scaleExtent([-20, 20])
    .on("zoom", () => {
      let e = d3.event.transform;
      e.x = currX;
      e.y = currY;
      if (e.k < currScale) {
        // center the svg when zoom-in
        currX = width / 2;
        currY = height / 2;
      }
      currScale = e.k;
      svg.attr(
        "transform",
        "translate(" + currX + "," + currY + ")scale(" + currScale + ")"
      );
    });
  var dragY, dragX;
  function dragstart() {
    let e = d3.event;
    dragX = e.x;
    dragY = e.y;
  }
  function dragend() {
    let e = d3.event;
    currX += e.x - dragX;
    currY += e.y - dragY;
    // let path = svg.select('.svgon')
    // let offset = 10
    // if ((Math.abs(e.x - dragX) < offset) && (Math.abs(e.y - dragY) < offset))
    //   return path.on('click').call(path.node(), path.datum())
  }
  function dragged() {
    let e = d3.event;
    let offX = currX + (e.x - dragX);
    let offY = currY + (e.y - dragY);
    svg.attr(
      "transform",
      "translate(" + offX + "," + offY + ")scale(" + currScale + ")"
    );
  }
  function mouseover(d) {
    // console.log('mouseover: ' + d.data.name + ': ' +
    //   d.centroid + ', angle: ' + d.textAngle + ', shown:' + d.shown)
    d3.select(this).classed("svgon", true);
    if (show_tooltips) {
      tooltip
        .transition()
        .duration(tooltipShowTime)
        .style("opacity", 0.9);
      tooltip
        .html(d.data.name + "<br/>" + formatNumber(d.value))
        .style("left", d3.event.clientX + "px")
        .style("top", d3.event.clientY - 70 + "px");
    }
  }
  function mouseout(d) {
    d3.select(this).classed("svgon", false);
    if (show_tooltips)
      tooltip
        .transition()
        .duration(tooltipHideTime)
        .style("opacity", 0);
  }
  function json2hierarchy(json) {
    let jh = d3.hierarchy(json).sort();
    jh.sum(d => {
      return d.root ? 0 : 1;
    });
    let partition = d3.partition();
    return partition(jh);
  }
  function reset_ui() {
    rootNode = clickedNode = null;
  }
  // d3.json("sunburst.json", function(err, json) {
  // if (err) throw err
  
  categories_json = json;
  set_hierarchy(json2hierarchy(json));
  reset_ui();
  d3.select(self.frameElement).style("height", height + "px");
  svg.call(zoom);
  // })

  var rootNode, clickedNode;
  function tween_labels(tr, clickArc) {
    if (!tr) tr = svg.transition().duration(tweenTransitionTime);
    tr
      .selectAll("text")
      .attrTween("transform", d => {
        return () => {
          return arcText(d);
        };
      })
      .attr("opacity", d => d.shown)
      .attr("height", 10)
      .attr("font-family", "monospace")
      .attr("font-size", function(d) {
        let id = d3
          .select(this)
          .attr("id")
          .substr(4);
        let bbox = d3
          .select("#path-" + id)
          .node()
          .getBBox();
        let bs = Math.min(bbox.width, bbox.height);
        let fs = Math.min(bs, d.outerR - d.innerR) / 4;
        let tl = this.getComputedTextLength() * 2;
        if (tl / fs > 10) fs /= 4;
        else if (tl / fs > 6) fs /= 2;
        return fs;
      })
      // .attr('text-anchor', (d) => (d.textAngle > 180 ? "start" : "end"))
      .attr("text-anchor", "middle");
  }
  function arcText(d) {
    var angle = (x((d.x0 + d.x1) / 2) - Math.PI / 2) / Math.PI * 180;
    d.textAngle = angle > 90 ? 180 + angle : angle;
    d.centroid = arc.centroid(d);
    if (!d.parent && !rootNode) rootNode = d;
    if (d == rootNode) d.centroid[0] = d.centroid[1] = d.textAngle = 0;
    else if (d.parent == rootNode && !d.addStrokeSpace) d.textAngle = 0;
    return "translate(" + d.centroid + ")rotate(" + d.textAngle + ")";
  }
  function click_tween(clickArc) {
    if (0) {
      // center-to-click
      currX += d3.event.clientX - currX;
      currY += d3.event.clientY - currY;
      svg.attr(
        "transform",
        "translate(" + currX + "," + currY + ")scale(" + currScale + ")"
      );
    }
    if (clickArc.parent) rootNode = clickArc.parent;
    clickedNode = clickArc;
    // console.log('click: ' + clickArc.data.name + ': ' +
    //  clickArc.centroid + ', angle: ' + clickArc.textAngle + ', shown:' + clickArc.shown)
    let tr = svg.transition().duration(tweenTransitionTime);
    tr.on("start", () => svg.selectAll("text").attr("opacity", 0));
    tr.on("end", () => tween_labels(tr, clickArc));
    tr.tween("scale", () => {
      var xd = d3.interpolate(x.domain(), [clickArc.x0, clickArc.x1]),
        yd = d3.interpolate(y.domain(), [clickArc.y0, 1]),
        yr = d3.interpolate(y.range(), [clickArc.y0 ? 40 : 0, radius]);
      return t => {
        x.domain(xd(t));
        y.domain(yd(t)).range(yr(t));
      };
    });
    tr
      .selectAll("path")
      .attr("opacity", d => {
        if (d == rootNode || (d.x0 >= clickArc.x0 && d.x1 <= clickArc.x1))
          d.shown = 1;
        else d.shown = 0;
        return d.shown;
      })
      .attrTween("d", d => {
        return () => arc(d);
      });
  }
  function set_hierarchy(hierarchy) {
    svg.selectAll("path").remove();
    svg.selectAll("text").remove();
    var paths = svg
      .selectAll("g")
      .data(hierarchy.descendants())
      .enter();
    paths
      .append("path")
      .attr("d", arc)
      .style("fill", d => {
        let n = d.children ? d : d.parent;
        if (!n) n = d;
        return color(n.data.name);
      })
      .on("mouseover", mouseover)
      .on("mouseout", mouseout)
      .on("click", click)
      .attr("id", (d, i) => "path-" + i);
    paths
      .append("text")
      .text(d => d.data.name)
      .attr("id", (d, i) => "tid-" + i);
    paths.exit().remove();
    // draw ui
    tween_labels();
  }
}
