function setipt(id, value) {
  Shiny.onInputChange(id, value)
}

$(document).ready(function() {
  Shiny.addCustomMessageHandler('background-color', function(color) {
    document.body.style.backgroundColor = color;
    //document.body.innerText = color;
  });
  Shiny.addCustomMessageHandler('computed', function(mess) {
    // The received value (in mess) is serialized in JSON, 
    // so we can  access the list element with object.name
    alert("Computed " + mess.what + " in " + mess.sec + " secs");
  })

  document.getElementById("mydiv").onclick = function() {
    var number = Math.random();
    Shiny.onInputChange("mydata", number);
    // document.getElementById("id").value = "value";
  };
  Shiny.addCustomMessageHandler("myCallbackHandler",
    function(color) {
      document.getElementById("mydiv").style.backgroundColor = color;
    }
  );
});

// Place title on the right of dashboard header. 
// function header() {
// $(document).ready(function() {
//   $("header").find("nav").append(\'<span class="myClass">spatialHeatmap</span>\');
// })
// }

// Full screen
function fullsn(elem) {
  if (elem.requestFullscreen) {
    elem.requestFullscreen();
  } else if (elem.mozRequestFullScreen) { /* Firefox */
    elem.mozRequestFullScreen();
  } else if (elem.webkitRequestFullscreen) { /* Chrome, Safari and Opera */
    elem.webkitRequestFullscreen();
  } else if (elem.msRequestFullscreen) { /* IE/Edge */
    elem.msRequestFullscreen();
  }
}

