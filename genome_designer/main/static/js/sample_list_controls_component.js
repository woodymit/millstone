/**
 * @fileoverview Component that decorates the controls for the list of
 *     ExperimentSamples.
 */


gd.SamplesControlsComponent = Backbone.View.extend({
  initialize: function() {
    this.listenToControls();

    this.maybeCreateFineS3Uploader();
  },

  listenToControls: function() {
    $('#gd-samples-create-from-server-location-submit').click(
        _.bind(this.handleCreateFromServerLocation, this));
  },

  handleCreateFromServerLocation: function() {
    // We validate by parsing the input elements, but then use the HTML5
    // FormData API to actually send the data to the server. FormData
    // makes the file upload "just work".
    var formDataForValidation = this.prepareRequestData(
        'gd-samples-create-from-server-location');
    if (!this.validateCreateFromServerLocation(formDataForValidation)) {
      return;
    }

    var onSuccess = function(responseData) {
      if (responseData.error.length) {
        alert('Error creating samples: ' + responseData.error);
        return;
      }

      // Success, reload the page.
      window.location.reload();
    }

    var formData = new FormData(
        $('#gd-samples-create-from-server-location')[0]);
    $.ajax({
      url: '/_/samples/create_from_server_location',
      type: 'POST',
      data: formData,
      success: onSuccess,

      // The following 3 param settings are necessary for properly passing
      // formData. See: http://stackoverflow.com/questions/166221/how-can-i-upload-files-asynchronously-with-jquery
      cache: false,
      contentType: false,
      processData: false
    });
  },

  /** Parses the form files and prepares the data. */
  prepareRequestData: function(formId) {
    var requestData = {}
    var formInputs = $('#' + formId + ' :input');
    _.each(formInputs, function(inputObj) {
      requestData[inputObj.name] = inputObj.value;
    });
    return requestData;
  },

  validateCreateFromServerLocation: function(requestData) {
    if (!requestData['targetsFile'].length) {
      alert('Please select a targest file to upload.');
      return false;
    }

    return true;
  },

  /**
   * Creates the S3 uploader if DOM target is present.
   *
   * TODO: Make sure this code still works.
   */
  maybeCreateFineS3Uploader: function() {
    if (this.$("#uploadDiv")) {
      this.uploader = this.$("#uploadDiv").fineUploaderS3({
        debug: true,
        multiple: true,
        autoUpload: true,
        maxConnections: 1,
        request: {
          endpoint: this.$("#uploadDiv").data("endpoint"),
          accessKey: this.$("#uploadDiv").data("accesskey")
        },
        signature: {
          endpoint: this.$("#uploadDiv").data("signature")
        },
        uploadSuccess: {
          endpoint: this.$("#uploadDiv").data("success")
        },
        objectProperties: {
          key: $.proxy(function(fileId) {
            var filename = this.uploader.fineUploader("getName", fileId);
            return "uploads/" + this.uuid + "/" + filename;
          }, this)
        },
        retry: {
          enableAuto: true
        },
        chunking: {
          enabled: true
        },
        deleteFile: {
          endpoint: this.$("#uploadDiv").data("delete"),
          enabled: true,
          forceConfirm: true
        },
        callbacks: {
          onError: function(id, name, reason) {
            alert(reason);
          }
        },
        validation: {
          stopOnFirstInvalidFile: false
        },
        template: 'samples-s3-uploader-template',
      }).on('complete', $.proxy(function(id, name, response, xhr) {
        var sid = xhr.s3file_id;

        // If we know sample_files, the newly uploaded file must be one of them. 
        if (this.sample_files) {
          this.sample_files[xhr.s3file_name]['sid'] = sid;

          not_uploaded = _.filter(_.values(this.sample_files), function(val) {
            return val.sid == undefined;
          })
          console.log("Missing " + not_uploaded.length + " sample files.");
          if (not_uploaded.length == 0) {
            $.post(this.$("#uploadDiv").data("finalize"), 
                JSON.stringify({
                  'sample_files': this.sample_files,
                  'targets_file_rows': this.targets_file_rows
                }),
                $.proxy(function(data) {
                  alert(data);
                }, this)
            );
          }
        }
        // The uploaded file is a template targets file that needs to be parsed by server.
        else {
          $.post(this.$("#uploadDiv").data("parse"), 
                {
                  's3file_id': sid,
                },
                $.proxy(function(data) {
                  console.log(data);
                  if (_.has(data, "error")) {
                    alert(data.error);
                  }
                  else {
                    // Save the sample files that we need from user as a dictionary
                    // and formatted rows of targets file to be sent along with s3file_ids
                    // when all samples are uploaded. 
                    this.targets_file_rows = data.targets_file_rows
                    this.sample_files = {}
                    _.each(data.sample_filenames, $.proxy(function(element, index, list) {
                      this.sample_files[element] = {};
                    }),this);
                    this.$(".qq-choose-file").text("Choose sample files ...")
                  }
                }, this)
          );
          // this.$("li.qq-upload-success").remove();
        }
      }, this)).on('validate', $.proxy(function(e, fileOrBlobData, button){
        // If we know sample_files, we must validate if dropped file exists in sample_files.
        if (this.sample_files) {
            if (_.has(this.sample_files, fileOrBlobData['name'])) {
              this.sample_files[fileOrBlobData['name']] = {'name': fileOrBlobData.name};
              return true;
            }
            else {
              console.log(fileOrBlobData['name'] + " is not listed in template targets file.")
              return false;
            }
        }
        else {
          return true;
        }
      }, this)).on('validateBatch', $.proxy(function(e, fileOrBlobDataArray, button){
        if (this.sample_files == undefined && fileOrBlobDataArray.length != 1) {
          alert("Cannot choose more than one template. ")
          return false;
        }
        return true;
      }, this));
    }
  },
});