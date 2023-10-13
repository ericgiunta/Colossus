#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstatsTODO {G1.0} *Statistical Software should list at least one primary reference from published academic literature.* 
#' @srrstatsTODO {G1.1} *Statistical Software should document whether the algorithm(s) it implements are:* - *The first implementation of a novel algorithm*; or - *The first implementation within **R** of an algorithm which has previously been implemented in other languages or contexts*; or - *An improvement on other implementations of similar algorithms in **R***. 
#' @srrstatsTODO {G1.3} *All statistical terminology should be clarified and unambiguously defined.* 
#' @srrstatsTODO {G1.5} *Software should include all code necessary to reproduce results which form the basis of performance claims made in associated publications.* 
#' @srrstatsTODO {G1.6} *Software should include code necessary to compare performance claims with alternative implementations in other R packages.* 
#' @srrstatsTODO {G2.0} *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
#' @srrstatsTODO {G2.0a} Provide explicit secondary documentation of any expectations on lengths of inputs
#' @srrstatsTODO {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#' @srrstatsTODO {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
#' @srrstatsTODO {G2.2} *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
#' @srrstatsTODO {G2.3} *For univariate character input:*
#' @srrstatsTODO {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @srrstatsTODO {G2.3b} *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
#' @srrstatsTODO {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstatsTODO {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstatsTODO {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#' @srrstatsTODO {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
#' @srrstatsTODO {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstatsTODO {G2.4e} *explicit conversion from factor via `as...()` functions*
#' @srrstatsTODO {G2.5} *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.* 
#' @srrstatsTODO {G2.6} *Software which accepts one-dimensional input should ensure values are appropriately pre-processed regardless of class structures.* 
#' @srrstatsTODO {G2.7} *Software should accept as input as many of the above standard tabular forms as possible, including extension to domain-specific forms.* 
#' @srrstatsTODO {G2.8} *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
#' @srrstatsTODO {G2.9} *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).* 
#' @srrstatsTODO {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.* 
#' @srrstatsTODO {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
#' @srrstatsTODO {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.* 
#' @srrstatsTODO {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' @srrstatsTODO {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.* 
#' @srrstatsTODO {G4.0} *Statistical Software which enables outputs to be written to local files should parse parameters specifying file names to ensure appropriate file suffices are automatically generated where not provided.* 
#' @srrstatsTODO {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*
#' @srrstatsTODO {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*
#' @srrstatsTODO {RE1.0} *Regression Software should enable models to be specified via a formula interface, unless reasons for not doing so are explicitly documented.*
#' @srrstatsTODO {RE1.1} *Regression Software should document how formula interfaces are converted to matrix representations of input data.* 
#' @srrstatsTODO {RE1.2} *Regression Software should document expected format (types or classes) for inputting predictor variables, including descriptions of types or classes which are not accepted.* 
#' @srrstatsTODO {RE1.3a} *Where otherwise relevant information is not transferred, this should be explicitly documented.* 
#' @srrstatsTODO {RE2.0} *Regression Software should document any transformations applied to input data, for example conversion of label-values to `factor`, and should provide ways to explicitly avoid any default transformations (with error or warning conditions where appropriate).*
#' @srrstatsTODO {RE2.1} *Regression Software should implement explicit parameters controlling the processing of missing values, ideally distinguishing `NA` or `NaN` values from `Inf` values (for example, through use of `na.omit()` and related functions from the `stats` package).* 
#' @srrstatsTODO {RE2.2} *Regression Software should provide different options for processing missing values in predictor and response data. For example, it should be possible to fit a model with no missing predictor data in order to generate values for all associated response points, even where submitted response values may be missing.*
#' @srrstatsTODO {RE3.0} *Issue appropriate warnings or other diagnostic messages for models which fail to converge.*
#' @srrstatsTODO {RE3.1} *Enable such messages to be optionally suppressed, yet should ensure that the resultant model object nevertheless includes sufficient data to identify lack of convergence.*
#' @srrstatsTODO {RE3.2} *Ensure that convergence thresholds have sensible default values, demonstrated through explicit documentation.*
#' @srrstatsTODO {RE3.3} *Allow explicit setting of convergence thresholds, unless reasons against doing so are explicitly documented.* 
#' @srrstatsTODO {RE4.0} *Regression Software should return some form of "model" object, generally through using or modifying existing class structures for model objects (such as `lm`, `glm`, or model objects from other packages), or creating a new class of model objects.*
#' @srrstatsTODO {RE4.1} *Regression Software may enable an ability to generate a model object without actually fitting values. This may be useful for controlling batch processing of computationally intensive fitting algorithms.* 
#' @srrstatsTODO {RE4.2} *Model coefficients (via `coef()` / `coefficients()`)*
#' @srrstatsTODO {RE4.3} *Confidence intervals on those coefficients (via `confint()`)*
#' @srrstatsTODO {RE4.4} *The specification of the model, generally as a formula (via `formula()`)*
#' @srrstatsTODO {RE4.5} *Numbers of observations submitted to model (via `nobs()`)*
#' @srrstatsTODO {RE4.6} *The variance-covariance matrix of the model parameters (via `vcov()`)*
#' @srrstatsTODO {RE4.7} *Where appropriate, convergence statistics* 
#' @srrstatsTODO {RE4.8} *Response variables, and associated "metadata" where applicable.*
#' @srrstatsTODO {RE4.9} *Modelled values of response variables.*
#' @srrstatsTODO {RE4.10} *Model Residuals, including sufficient documentation to enable interpretation of residuals, and to enable users to submit residuals to their own tests.*
#' @srrstatsTODO {RE4.11} *Goodness-of-fit and other statistics associated such as effect sizes with model coefficients.*
#' @srrstatsTODO {RE4.12} *Where appropriate, functions used to transform input data, and associated inverse transform functions.* 
#' @srrstatsTODO {RE4.13} *Predictor variables, and associated "metadata" where applicable.* 
#' @srrstatsTODO {RE4.14} *Where possible, values should also be provided for extrapolation or forecast *errors*.*
#' @srrstatsTODO {RE4.15} *Sufficient documentation and/or testing should be provided to demonstrate that forecast errors, confidence intervals, or equivalent values increase with forecast horizons.* 
#' @srrstatsTODO {RE4.16} *Regression Software which models distinct responses for different categorical groups should include the ability to submit new groups to `predict()` methods.* 
#' @srrstatsTODO {RE4.17} *Model objects returned by Regression Software should implement or appropriately extend a default `print` method which provides an on-screen summary of model (input) parameters and (output) coefficients.*
#' @srrstatsTODO {RE4.18} *Regression Software may also implement `summary` methods for model objects, and in particular should implement distinct `summary` methods for any cases in which calculation of summary statistics is computationally non-trivial (for example, for bootstrapped estimates of confidence intervals).* 
#' @srrstatsTODO {RE5.0} *Scaling relationships between sizes of input data (numbers of observations, with potential extension to numbers of variables/columns) and speed of algorithm.* 
#' @srrstatsTODO {RE6.0} *Model objects returned by Regression Software (see* **RE4***) should have default `plot` methods, either through explicit implementation, extension of methods for existing model objects, or through ensuring default methods work appropriately.*
#' @srrstatsTODO {RE6.1} *Where the default `plot` method is **NOT** a generic `plot` method dispatched on the class of return objects (that is, through an S3-type `plot.<myclass>` function or equivalent), that method dispatch (or equivalent) should nevertheless exist in order to explicitly direct users to the appropriate function.*
#' @srrstatsTODO {RE6.2} *The default `plot` method should produce a plot of the `fitted` values of the model, with optional visualisation of confidence intervals or equivalent.* 
#' @srrstatsTODO {RE6.3} *Where a model object is used to generate a forecast (for example, through a `predict()` method), the default `plot` method should provide clear visual distinction between modelled (interpolated) and forecast (extrapolated) values.* 
#' @srrstatsTODO {RE7.3} Demonstrate and test expected behaviour when objects returned from regression software are submitted to the accessor methods of **RE4.2**--**RE4.7**.
#' @srrstatsTODO {RE7.4} Extending directly from **RE4.15**, where appropriate, tests should demonstrate and confirm that forecast errors, confidence intervals, or equivalent values increase with forecast horizons.
#' @noRd
NULL

#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#' @srrstatsNA {G3.1} Only one method avaliable
#' @srrstatsNA {G3.1a} Arbitrary methods not explicitly supported
#' @srrstatsNA {G3.0} Event times must be compared, order is significant so a tolerance would not be appropriate unless the user performs rounding 
#' @srrstatsNA {G5.9} Some noise CAN change results
#' @srrstatsNA {G5.9a} Noise in data is not always trivial
#' @srrstatsNA {G5.10} No extended tests used
#' @srrstatsNA {G5.11} Size of examples used in publication will be provided there
#' @srrstatsNA {G5.11a} No data downloaded
#' @srrstatsNA {G5.12} Conditions will be provided when data is provided, none known for current testing
#' @srrstatsNA {RE1.3} Most functions correct the order of options, missing control parameters, or output values with no aspects to retain
#' @srrstatsNA {RE2.3} Centering not relevent
#' @srrstatsNA {RE2.4} Collinearity not intended to be detected
#' @srrstatsNA {RE2.4a} *Perfect collinearity among predictor variables*
#' @srrstatsNA {RE2.4b} Perfect collinearity between covariates and events exceptionally improbable
#' @srrstatsNA {RE7.0} All data has noise
#' @srrstatsNA {RE7.0a} noiseless data isn't refused
#' @srrstatsNA {RE7.1} All data has noise
#' @srrstatsNA {RE7.1a} Not relevent to compare noisy to non-noisy data
#' @srrstatsNA {RE7.2} Tests verify that inputs are only changed as intended
#' @noRd
NULL
