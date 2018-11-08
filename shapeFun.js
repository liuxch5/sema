

shinyjs.printGraph = function(){
    svg = new yfiles.view.SvgExport(graphControl.viewport, 2);
    element = svg.exportSvg(graphControl);
    html = yfiles.view.SvgExport.exportSvgString(element);
    const w = window.open();
   
    w.document.write(html);
    w.document.close();
    w.print();
    w.close();
};

shinyjs.toImage = function(){
    svg = new yfiles.view.SvgExport(graphControl.viewport, 2);
    svgElement = svg.exportSvg(graphControl);
    svgElement.style.overflow = "scroll";
    html_svg = yfiles.view.SvgExport.exportSvgString(svgElement);
    Shiny.onInputChange("jsImage", html_svg);

};

shinyjs.clearEdges = function () {
    if(graph.edges !== null) {
        var remEdg = [];
        for(var i = 0; i < graph.edges.size; i++){
            remEdg[i] = graph.edges.get(i);
            // graph.remove(graph.edges.get(i));

        }

        for (var i = 0; i < remEdg.length; i++) {
            graph.remove(remEdg[i]);
        }
    }
};

var lay;
shinyjs.layout = function (params) {
    //var tmp = typeof params;
    // console.log(params);
    if(params === null)
        shinyjs.layout("Tree");
    if(graph!==null){
        if((typeof params) !== "String"){
            params = params.toString();
        }
        //alert("new " + typeof params)
        switch (params){
            case "Hierarchical":
                lay = new yfiles.hierarchic.HierarchicLayout();
                lay.edgeLayoutDescriptor.minimumDistance = 80;
                graphControl.morphLayout(lay);
                break;
            case "Organic":
                lay = new yfiles.organic.OrganicLayout();
                lay.preferredEdgeLength = 150;
                graphControl.morphLayout(lay);
                break;
            case "Orthogonal":
                graphControl.morphLayout(new yfiles.orthogonal.OrthogonalLayout());
                break;
            case "Tree":
                graphControl.morphLayout(new yfiles.tree.TreeLayout());
                break;
            case "Radial":
                graphControl.morphLayout(new yfiles.radial.RadialLayout());
                break;
        }
        graphControl.fitGraphBounds();
    }



};

shinyjs.submit = function () {
    var nodes = [];
    var edges = [];
    var eStyles = [];
    var types = [];
    var node;
    var style;
    var edge;

    var nodeArray = graph.nodes;
    var edgeArray = graph.edges;

    for (i = 0; i < nodeArray.size; i++) {
        node = nodeArray.elementAt(i);
        style = node.tag;
        nodes[i] = style.label;
        types[i] = style.type;
    }

    for (i = 0; i < edgeArray.size; i++) {
        edge = edgeArray.elementAt(i);
        edges[i] = edge.sourceNode.tag.label + ":" + edge.targetNode.tag.label;
        eStyles[i] = edge.tag;
    }
    
    var tmp;
    for (var i = 0; i < eStyles.length; i++) {
        if (eStyles[i] === "Factor") {
            tmp = edges[i].split(":");
            tmp = tmp[1].split(".");
        }
    }


    var g = {
        nodes: nodes,
        edges: edges,
        edgeTypes: eStyles
    };
    Shiny.onInputChange("graph", g);
};

shinyjs.edgeChange = function (params) {
    edgeType = params[0];
    var style = new yfiles.styles.ArcEdgeStyle();
    style.sourceArrow = yfiles.styles.IArrow.NONE;
    if (edgeType === "Causal")
        style.targetArrow = yfiles.styles.IArrow.DEFAULT;
    else
        style.targetArrow = yfiles.styles.IArrow.CIRCLE;
    style.stroke = new yfiles.view.Stroke(new yfiles.view.SolidColorFill(new yfiles.view.Color(0, 0, 0, 255)), 1);
    graph.edgeDefaults.style = style;
};

shinyjs.results = function (params) {
    // check if the flag is up
    // if so clear the graph
    if(params === true){
        console.log("remove results");
        shinyjs.clearEdges();
        return;
    }
       
    var edgeArray = graph.edges;
    console.log(edgeArray.size + "\n");
    lhs = params[0].lhs;
    rhs = params[0].rhs;
    std = params[0].std;
    val = params[0].val;

    for (i = 0; i < edgeArray.size; i++) {
        edge = edgeArray.elementAt(i);
        node1 = edge.targetNode.tag.label;
        node2 = edge.sourceNode.tag.label;
        console.log(node1 + " " + node2);

        for (j = 0; j < lhs.length; j++) {
            if (lhs[j] === node1 && rhs[j] === node2) {
                col = yfiles.view.Fill.BLACK;
                if (std[j] > 0)
                    col = yfiles.view.Fill.RED;
                else if (std[j] < 0)
                    col = yfiles.view.Fill.BLUE;
                style = new yfiles.styles.ArcEdgeStyle();
                style.sourceArrow = edge.style.sourceArrow;
                style.targetArrow = edge.style.targetArrow;
                style.stroke = new yfiles.view.Stroke(col, 0.5 + Math.abs(std[j]) * 10);
                graph.setStyle(edge, style);
                if(val[j] != "") graph.setLabelText(edge.labels.get(0), parseFloat(val[j]).toFixed(3));
                                
                break;
            }
        }
    }    
    graphControl.updateVisual();    
};



var flag1 = true, flag2 = false;

shinyjs.addNode = function (params) {

    //window.alert(graph.nodes);
    
    //window.alert(params.label);

    if (params.label === "") {
        if (flag1 === true)
            flag1 = false;
        else if(flag1 === false)
            flag1 = true;
        Shiny.onInputChange("flagJS2", flag1);
        //window.alert("Node name cannot be blank...");
        return;
    }

    var sns = new yfiles.styles.ShapeNodeStyle();

    var dp = Object.create(NodeStyle);

    pr = shinyjs.getParams(params, dp);
    

    if (pr.type === "mRNA (RNAseq)") {
        sns.shape = yfiles.styles.ShapeNodeShape.RECTANGLE;
        pr.label = pr.label + ".RNA";
    }
    if (pr.type === "Log-mRNA (RNAseq)") {
        sns.shape = yfiles.styles.ShapeNodeShape.RECTANGLE;
        pr.label = pr.label + ".LogRNA";
    }
    if (pr.type === "Protein (RPPA)") {
        sns.shape = yfiles.styles.ShapeNodeShape.DIAMOND;
        pr.label = pr.label + ".RPPA";
    }
    if (pr.type === "Copy Number Variation (SNP6)") {
        sns.shape = yfiles.styles.ShapeNodeShape.HEXAGON;
        pr.label = pr.label + ".CNV";
    }
    if (pr.type === "Somatic mutation (WXS)") {
        sns.shape = yfiles.styles.ShapeNodeShape.TRAPEZ2;
        pr.label = pr.label + ".Mut";
    }
    if (pr.type === "Factor") {
        sns.shape = yfiles.styles.ShapeNodeShape.ELLIPSE;
        pr.label = pr.label + ".Factor";
    }
    if (pr.type === "miRNA (miRNA-seq)") {
        sns.shape = yfiles.styles.ShapeNodeShape.TRIANGLE;
        pr.label = pr.label + ".miRNA";
    }
    if (pr.type === "Tumor features (PanCan iAtlas)") {
        sns.shape = yfiles.styles.ShapeNodeShape.OCTAGON;
        pr.label = pr.label + ".iAtlas";
    }
    if (pr.type === "Germline variation (WXS)") {
        sns.shape = yfiles.styles.ShapeNodeShape.STAR8;
        pr.label = pr.label + ".Germ";
    }
    if (pr.type === "Clinical") {
        sns.shape = yfiles.styles.ShapeNodeShape.STAR5;
        if(pr.label === "Overall_Survival") pr.label = pr.label + ".Surv";
        else pr.label = pr.label + ".Clin";
    }

    if (pr.type === "Group") {
        sns.shape = yfiles.styles.ShapeNodeShape.TRIANGLE2;        
    }

    for (var i = 0; i < graph.labels.size; i++) {

        if (pr.label === graph.labels.get(i).text) {
            if (flag2 === true)
                flag2 = false;
            else if(flag2 === false)
                flag2 = true;
            Shiny.onInputChange("flagJS1", flag2);
            return;
        }
    }    
    var node = graph.createNode(new yfiles.geometry.Rect(100, 200, 50, 50), sns, pr);
    graph.addLabel(node, pr.label);    
};

findNode = function(label){    
    for(var i = 0; i < graph.nodes.size; i++){
        lab = graph.nodes.get(i).tag.label;
        if(lab === label) return(graph.nodes.get(i));
    }
    return(null);
};

shinyjs.addEdge = function(params){
    
    var dp = Object.create(EdgeStyle);
    pr = shinyjs.getParams(params, dp);
    label1 = findNode(pr.left);
    label2 = findNode(pr.right);
    type = pr.type;

    var style = new yfiles.styles.ArcEdgeStyle();
    style.sourceArrow = yfiles.styles.IArrow.NONE;
    if (pr.type === "Causal")
        style.targetArrow = yfiles.styles.IArrow.DEFAULT;
    else
        style.targetArrow = yfiles.styles.IArrow.CIRCLE;
    style.stroke = new yfiles.view.Stroke(new yfiles.view.SolidColorFill(new yfiles.view.Color(0, 0, 0, 255)), 1);

    var edge = graph.createEdge(label1, label2, style);
};

shinyjs.copyGraph = function(){
    clone = new yfiles.graph.DefaultGraph();
    copier = new yfiles.graph.GraphCopier();
    copier.clone = yfiles.graph.CloneTypes.ALL;
    copier.copy(graph, clone);
    graphs[graphs.length] = clone;
    console.log("Cloning graph " + graphs.length);
};

shinyjs.loadGraph = function(params){
    console.log("Model " + params);
    clone = new yfiles.graph.DefaultGraph();
    copier = new yfiles.graph.GraphCopier();
    copier.clone = yfiles.graph.CloneTypes.ALL;
    copier.copy(graphs[parseInt(params) - 1], clone);    
    graph = clone;

    graphControl.graph = graph;
    var et = [document.getElementById("edgeSelect").value, "sd"];
    console.log(et);
    
    shinyjs.edgeChange(et);
    graph.addEdgeCreatedListener(function(object, event){

        event.item.tag = edgeType;
        layoutType = document.getElementById("layoutSelect").value;
        //alert(typeof layoutType);
        shinyjs.layout(layoutType);
        graph.addLabel(event.item, "");

    });
    graphControl.updateVisual();
};

shinyjs.clearGraph = function(){
    graph.clear();
};

shinyjs.renderLabel = function (params) {
    // If the params is empty don't show any label on edges
    lhs = params[0].lhs;
    rhs = params[0].rhs;
    val = params[0].val;
        
    //mark = params[0][" + mark + "]
    combTable = [];
    combEdges = [];

    for(var i =0; i < lhs.length; i++) {
        combTable[i] = lhs[i] + ":" + rhs[i];
    }

    for(var i = 0; i < graph.edges.size; i++){
        combEdges[i] = graph.edges.get(i).targetNode.tag.label + ":" + graph.edges.get(i).sourceNode.tag.label;
    }

    for(var i = 0; i < combTable.length; i++){
        for(var j = 0; j < combEdges.length; j++){
            if(combTable[i] === combEdges[j]){
                if(isNumber(val[i]))
                    graph.setLabelText(graph.edges.get(j).labels.get(0), parseFloat(val[i]).toFixed(3));
                else graph.setLabelText(graph.edges.get(j).labels.get(0), val[i]);
            }
        }
    }
    graphControl.updateVisual();
};

isNumber = function(n) {
  return !isNaN(parseFloat(n));
};

// focuses to text box for naming variable
shinyjs.focus = function(){
    var selectized = $('#dynamic').selectize();
    selectized[0].selectize.focus();
};

