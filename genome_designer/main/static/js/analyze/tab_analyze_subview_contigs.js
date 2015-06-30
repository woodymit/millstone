/**
 * @fileoverview Sub-view displaying Variants.
 */


gd.TabAnalyzeSubviewContigs = gd.TabAnalyzeSubviewAbstractBase.extend(
{
  /** @override */
  initialize: function() {
    this.render();
  },

  render: function() {
    this.redrawDatatable();
  },

  /** Draws or redraws the table. */
  redrawDatatable: function() {
    if (this.datatableComponent) {
      this.datatableComponent.destroy();
    }

    var requestData = {
        refGenomeUid: this.model.attributes.refGenomeUid,
        alignmentGroupUid: this.model.attributes.alignmentGroupUid
    };

    this.datatableComponent = new gd.DataTableComponent({
        el: $('#gd-datatable-hook'),
        serverTarget: '/_/contigs',
        controlsTemplate: '/_/templates/contig_list_controls',
        requestData: requestData,
    });

    this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
        _.bind(this.decorateControls, this));
  },

  decorateControls: function() {
    this.contigControlsComponent = new gd.ContigControlsComponent({
      el: '#gd-datatable-hook-control',
      datatableComponent: this.datatableComponent,
      alignmentGroupUid: this.model.attributes.alignmentGroupUid
    });

    this.listenTo(this.contigControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawDatatable, this));
  }
});
