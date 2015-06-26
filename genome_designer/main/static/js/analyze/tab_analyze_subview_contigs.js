/**
 * @fileoverview Sub-view displaying Variants.
 */


gd.TabAnalyzeSubviewContigs = gd.TabAnalyzeSubviewAbstractBase.extend(
{
  /** @override */
  initialize: function() {
    this.variantsTableComponent = null;

    // List of variant sets belonging to the user.
    this.variantSetList = null;

    // Map of variant keys, used to draw the fields dropdown.
    this.variantKeyMap = null;

    // Keys that should be visible.
    // When null, we are just showing the default.
    this.visibleKeyNames = null;

    // Check whether a filter is set on the view.
    var paramObject = gd.Util.getQueryStringAsObject();
    var maybeQueryString = window.location.search;
    if ('filter' in paramObject) {
      var filterString = paramObject['filter'];
      this.model.set('filterString', filterString);
    } else {
      this.model.set('filterString', '');
    }

    // Set the view to cast as default, unless melt=1 param in url.
    this.model.set('isMelted',
        'melt' in paramObject && paramObject['melt'] == '1');

    this.variantsTableComponent = new gd.VariantsTableComponent({
      model: this.model
    });

    this.setupListeners();
  },

  setupListeners: function() {
    // Propagate NAVIGATE event.
    // TODO: Is there more elegant way to do this?
    this.listenTo(this.variantsTableComponent, 'NAVIGATE', function (e) {
      this.trigger('NAVIGATE', e);
    });
  },

  // /** @override */
  // destroy: function() {
  //   if (this.variantsTableComponent) {
  //     this.variantsTableComponent.destroy()
  //   }
  //   this.variantsTableComponent = null;
  // }

  // /** Draws or redraws the table. */
  // redrawDatatable: function() {
  //   if (this.datatableComponent) {
  //     this.datatableComponent.destroy();
  //   }

  //   var requestData = {
  //       projectUid: this.model.get('projectUid'),
  //       showDeNovo: this.showDeNovo ? 1 : 0
  //   };

  //   this.datatableComponent = new gd.DataTableComponent({
  //       el: $('#gd-datatable-hook'),
  //       serverTarget: '/_/ref_genomes',
  //       controlsTemplate: '/_/templates/reference_genome_list_controls',
  //       requestData: requestData,
  //   });

  //   this.listenTo(this.datatableComponent, 'DONE_CONTROLS_REDRAW',
  //       _.bind(this.decorateControls, this));
  // },

  // decorateControls: function() {
  //   this.refGenomeControlsComponent = new gd.RefGenomeControlsComponent({
  //     el: '#gd-ref-genome-list-view-datatable-hook-control',
  //     datatableComponent: this.datatableComponent
  //   });

  //   this.listenTo(this.refGenomeControlsComponent, 'MODELS_UPDATED',
  //       _.bind(this.redrawDatatable, this));
  //   this.listenTo(this.refGenomeControlsComponent, 'TOGGLE_DE_NOVO',
  //       _.bind(function() {
  //           this.showDeNovo = !this.showDeNovo;
  //           this.redrawDatatable();
  //       }, this));
  // }
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

