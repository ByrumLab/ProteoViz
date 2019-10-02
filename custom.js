document.getElementById("EGSEAplot1button").addEventListener("click", myFunction3);

function myFunction3() {
  var plotHeight = document.getElementById("EGSEA1_plot_height").value + "px";
  document.getElementById("EGSEA1").style = "height:" + plotHeight;
}
