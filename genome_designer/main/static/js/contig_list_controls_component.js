/**
 * @fileoverview Component that decorates the controls for the list of
 *     ReferenceGenomes.
 */


gd.ContigControlsComponent = gd.DataTableControlsComponent.extend({
  initialize: function() {
    gd.DataTableControlsComponent.prototype.initialize.call(this);

    this.decorateControls();
  },

  decorateControls: function() {
    $('#gd-contig-assemble-submit').click(
        _.bind(this.handleAssembleContigs, this));

    this.drawDropdownOptions();
  },

  /** Puts UI in the loading state. */
  enterLoadingState: function() {
    $(".gd-id-form-submit-button")
        .prop('disabled', true);

    this.loadingSpinner = new gd.Spinner();
    this.loadingSpinner.spin();
  },

  /** Puts UI in the loading state. */
  exitLoadingState: function() {
    $(".gd-id-form-submit-button")
        .prop('disabled', false);
    this.loadingSpinner.stop();
  },

   /** Draws dropdown options. */
  drawDropdownOptions: function() {
    // Option to delete samples.
    var deleteOptionHtml =
        '<a href="#" class="gd-id-contigs-delete">Delete</a>';
    this.addDropdownOption(deleteOptionHtml);
    $('.gd-id-contigs-delete').click(_.bind(this.handleDelete, this));
  },

  // drawNewDropdownOptions: function() {
  //   // Generate contigs with default parameters
  //   $('.gd-id-generate-contigs-default').click(_.bind(this.handleAssembleContigs, this));
  // },

  /** Send request to generate contigs with default parameters **/
  handleAssembleContigs: function() {
    var requestData = this.prepareRequestData('gd-contig-assemble');
    console.log('requestData:')
    console.log(requestData)
  },

  /** Parses the form files and prepares the data. */
  prepareRequestData: function(formId) {
    var requestData = {}

    // This works correctly for all inputs except radio buttons
    var formInputs = $('#' + formId + ' :input');
    _.each(formInputs, function(inputObj) {
      requestData[inputObj.name] = inputObj.value;
    });

    return requestData;
  },

  /** Sends request to delete selected samples. */
  handleDelete: function() {
    var contigUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!contigUidList.length) {
      alert('Please select contigs to delete.');
      return;
    }

    // Get confirmation from user.
    var agree = confirm('Are you sure you want delete these contigs?')
    if (!agree) {
      return;
    }

    this.enterLoadingState();

    var postData = {
        contigUidList: contigUidList,
    };

    $.post('/_/contigs/delete', JSON.stringify(postData),
        _.bind(this.handleDeleteResponse, this));
  },

  handleDeleteResponse: function(response) {
    this.exitLoadingState();

    if ('error' in response && response.error.length) {
      alert(response.error);
    } else {
      this.trigger('MODELS_UPDATED');
    }
  },

  handleConcatenate: function() {
    var refGenomeUidList = this.datatableComponent.getCheckedRowUids();

    // If nothing to do, show message.
    if (!refGenomeUidList.length) {
      alert('Please select reference genomes to concatenate.');
      return;
    }
    // If only one selected, show message.
    if (refGenomeUidList.length == 1) {
      alert('Please select more than one reference genome to concatenate.');
      return;
    }

    // Get new genome name
    var newGenomeLabel = prompt(
        'Enter a name for the concatenated genome:', 'new_genome_name');
    while (newGenomeLabel == '') {
      var newGenomeLabel = prompt(
          'Please enter a non-zero length name for the concatenated genome',
          'new_genome_name');
    }
    if (newGenomeLabel == null) {
      return;
    }

    this.enterLoadingState();

    var postData = {
        'newGenomeLabel': newGenomeLabel,
        'refGenomeUidList': refGenomeUidList,
    };

    $.post('/_/ref_genomes/concatenate', {data:JSON.stringify(postData)},
        _.bind(this.handleConcatenateResponse, this));

},

  handleConcatenateResponse: function(response) {
    this.exitLoadingState();

    if ('error' in response && response.error.length) {
      alert(response.error);
    } else {
      this.trigger('MODELS_UPDATED');
    }
},
});
