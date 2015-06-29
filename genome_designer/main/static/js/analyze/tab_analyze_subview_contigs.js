/**
 * @fileoverview Sub-view displaying Variants.
 */


gd.TabAnalyzeSubviewContigs = gd.TabAnalyzeSubviewAbstractBase.extend(
{
  /** @override */
  initialize: function() {
    console.log('HELLO')
    this.showDeNovo = true;
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
    console.log('this.datatableComponent:')
    console.log(this.datatableComponent)

    this.listenTo(this.contigControlsComponent, 'MODELS_UPDATED',
        _.bind(this.redrawDatatable, this));
  }
});


// /**
//  * @fileoverview Reference Genome List view.
//  */


// gd.RefGenomeListView = Backbone.View.extend({
//   el: '#gd-page-container',

//   initialize: function() {
//     this.showDeNovo = true;
//     this.render();
//   },

//   render: function() {
//     $('#gd-sidenav-link-refgenomes').addClass('active');

//     this.redrawDatatable();
//   },

//   /** Draws or redraws the table. */
//   redrawDatatable: function() {
//     if (this.datatableComponent) {
//       this.datatableComponent.destroy();
//     }

//     var requestData = {
//         projectUid: this.model.get('uid'),
//         showDeNovo: this.showDeNovo ? 1 : 0
//     };

//     this.datatableComponent = new gd.DataTableComponent({
//         el: $('#gd-datatable-hook'),
//         serverTarget: '/_/ref_genomes',
//         controlsTemplate: '/_/templates/reference_genome_list_controls',
//         requestData: requestData,
//     });

//     this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
//         _.bind(this.decorateControls, this));
//   },

//   decorateControls: function() {
//     this.refGenomeControlsComponent = new gd.RefGenomeControlsComponent({
//       el: '#gd-ref-genome-list-view-datatable-hook-control',
//       datatableComponent: this.datatableComponent
//     });

//     this.listenTo(this.refGenomeControlsComponent, 'MODELS_UPDATED',
//         _.bind(this.redrawDatatable, this));
//     this.listenTo(this.refGenomeControlsComponent, 'TOGGLE_DE_NOVO',
//         _.bind(function() {
//             this.showDeNovo = !this.showDeNovo;
//             this.redrawDatatable();
//         }, this));
//   }
// });

