/**
 * Created by KOM7NL on 12/2/2016.
 */

var graphControl;
var graph;
var edge;
var graphs = [];

let fileSupport = null;

var NodeStyle = {
    label: "",
    type: "",
    index: 0
};
var EdgeStyle = {
    left: "",
    right: "",
    type: ""    
};
var resStyle = {
    lhs: [],
    rhs: [],
    z: [],
    std: []
};
var edgeType = "Causal";
var size;
var target;
var pass;

graphControl = new yfiles.view.GraphComponent("#graphComponent");
graphControl.mouseWheelBehavior = yfiles.view.MouseWheelBehaviors.ZOOM;    

const igraph = graphControl.graph;
graph = igraph.lookup(yfiles.graph.DefaultGraph.$class);
//console.log(graph);
var mode = new yfiles.input.GraphEditorInputMode();

//mode.allowCreateBend = false;
mode.allowCreateNode = false;
mode.allowEditLabel = false;
graphControl.inputMode = mode;
graphControl.graph.nodeDefaults.size = new yfiles.geometry.Size(50, 50);
var sns = new yfiles.styles.ShapeNodeStyle();
sns.shape = yfiles.styles.ShapeNodeShape.RECTANGLE;
graph.nodeDefaults.style = sns;
var es = new yfiles.styles.ArcEdgeStyle();
es.stroke = new yfiles.view.Stroke(new yfiles.view.SolidColorFill(new yfiles.view.Color(0, 0, 0, 255)), 1);
es.sourceArrow = yfiles.styles.IArrow.NONE;
es.targetArrow = yfiles.styles.IArrow.DEFAULT;    

graph.edgeDefaults.style = es;
var createEdge = mode.createEdgeInputMode;
var layoutType;
graph.addEdgeCreatedListener((object, event) => {

    event.item.tag = edgeType;
    layoutType = document.getElementById("layoutSelect").value;
    //alert(typeof layoutType);
    shinyjs.layout(layoutType);
    graph.addLabel(event.item, "");

});

createEdge.addEdgeCreatedListener((object, event) => {
    target = event.item.targetNode.toString().split(".");
    source = event.item.sourceNode.toString().split(".");
    if(source[1] == "Surv"){
        graph.remove(event.item);
    } else if(edgeType == "Factor" & target[1] != "Factor"){            
        graph.remove(event.item);            
    } 
});

mode.addDeletedItemListener((object, event) => {
    if(event.item.tag != null){
        size = Object.keys(event.item).length;

        // If a node has no edge object size is 9; if it has edge/edges the size is 12
        if (size == 12 || size == 9) {            
            console.log(event.item.tag.label);
            Shiny.onInputChange("nodeDeleted", event.item.tag.label);
        }        
    }     
});

mode.addItemClickedListener((object, event) => {
    // if item.tag is null edge label is clicked
    if(event.item.tag != null){
        size = Object.keys(event.item).length;

        // If a node has no edge object size is 9; if it has edge/edges the size is 12
        if (size == 12 || size == 9) {
            pass = String(event.item);
            Shiny.onInputChange("clicked", pass);
        }
        //Edge size is 15
        else if (size == 15) {
            //window.alert(typeof event)
            // graph.setLabelText(event.item, "hello");
            pass = String(event.item.sourceNode) + " : " + String(event.item.targetNode);
            //alert(pass)
            Shiny.onInputChange("clicked", pass);
        }
    }
});
    
mode.addMultiSelectionFinishedListener((object, event) =>{
    selection = graphControl.selection.selectedNodes.toArray();
    res = String(selection);
    Shiny.onInputChange("multiSelection", res);
});


    


