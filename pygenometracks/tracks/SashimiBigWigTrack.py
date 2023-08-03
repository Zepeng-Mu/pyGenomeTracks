from .GenomeTrack import GenomeTrack
import numpy as np
from ..utilities import plot_coverage, InputError, transform, change_chrom_names, opener, to_string, change_chrom_names, temp_file_from_intersect, get_region
import pyBigWig
from intervaltree import IntervalTree, Interval
import matplotlib
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from tqdm import tqdm

Path = mpath.Path

DEFAULT_LINKS_COLOR = 'blue'
HUGE_NUMBER = int(1e9)  # Which should be above any chromosome size

DEFAULT_BIGWIG_COLOR = '#33a02c'


class SashimiBigWigTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.bw', '.bigwig', '.bam', '.sashimi']
    TRACK_TYPE = 'sashimiBigWig'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
color = #666666
# To use a different color for negative values
#negative_color = red
# To use transparency, you can use alpha
# default is 1
# alpha = 0.5
# the default for min_value and max_value is 'auto' which means that the scale will go
# roughly from the minimum value found in the region plotted to the maximum value found.
min_value = 0
#max_value = auto
# The number of bins takes the region to be plotted and divides it
# into the number of bins specified
# Then, at each bin the bigwig mean value is computed and plotted.
# A lower number of bins produces a coarser tracks
number_of_bins = 700
# to convert missing data (NaNs) into zeros. Otherwise, missing data is not plotted.
nans_to_zeros = true
# The possible summary methods are given by pyBigWig:
# mean/average/stdev/dev/max/min/cov/coverage/sum
# default is mean
summary_method = mean
# for type, the options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similarly points:ms sets the point size (markersize (ms) to the given float
# type = line:0.5
# type = points:0.5
# set show_data_range to false to hide the text on the left showing the data range
show_data_range = true
# to compute operations on the fly on the file
# or between 2 bigwig files
# operation will be evaluated, it should contains file or
# file and second_file,
# we advice to use nans_to_zeros = true to avoid unexpected nan values
#operation = 0.89 * file
#operation = - file
#operation = file - second_file
#operation = log2((1 + file) / (1 + second_file))
#operation = max(file, second_file)
#second_file = path for the second file
# To log transform your data you can also use transform and log_pseudocount:
# For the transform values:
# 'log1p': transformed_values = log(1 + initial_values)
# 'log': transformed_values = log(log_pseudocount + initial_values)
# 'log2': transformed_values = log2(log_pseudocount + initial_values)
# 'log10': transformed_values = log10(log_pseudocount + initial_values)
# '-log': transformed_values = - log(log_pseudocount + initial_values)
# For example:
#tranform = log
#log_pseudocount = 2
# When a transformation is applied, by default the y axis
# gives the transformed values, if you prefer to see
# the original values:
#y_axis_values = original
# If you want to have a grid on the y-axis
#grid = true
# The link in Sashimi plot is a Bezier curve.
# The height of the curve is calculated from the length of the intron.
# When the y-axis in bigwig track is different, the height of curve needs to be scaled.
scale_link_height = 1
# The line width for links is proportion to the numbers at the last column in links file (PSI).
# But the absolute width is calculated from the supplied numbers, which can look too thin or too wide sometimes.
# Use scale_line_width to scale the absolute line widths.
# You may need to try several values to get a satisfying result.
# Use this to tell pyGenomeTracks whether to label PSI on links
show_number = false
scale_line_width = 3
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = {
        'max_value': None,
        'min_value': None,
        'show_data_range': True,
        'orientation': None,
        'bw_color': DEFAULT_BIGWIG_COLOR,
        'link_color': DEFAULT_LINKS_COLOR,
        'negative_color': None,
        'alpha': 1,
        'nans_to_zeros': False,
        'summary_method': 'mean',
        'number_of_bins': 700,
        'type': 'fill',
        'transform': 'no',
        'log_pseudocount': 0,
        'y_axis_values': 'transformed',
        'second_file': None,
        'operation': 'file',
        'grid': False,
        'line_width': None,
        'line_style': 'solid',
        'max_value': None,
        'min_value': None,
        'show_number': False,
        'region': None,  # Cannot be set manually but is set by tracksClass
        'ylim': None,
        'use_middle': False,
        'scale_link_height': 1,
        'scale_line_width': 1,
    }
    NECESSARY_PROPERTIES = ['file', 'link_file']
    SYNONYMOUS_PROPERTIES = {
        'max_value': {
            'auto': None
        },
        'min_value': {
            'auto': None
        },
    }
    POSSIBLE_PROPERTIES = {
        'orientation': [None, 'inverted'],
        'summary_method': [
            'mean', 'average', 'max', 'min', 'stdev', 'dev', 'coverage', 'cov',
            'sum'
        ],
        'transform': ['no', 'log', 'log1p', '-log', 'log2', 'log10'],
        'y_axis_values': ['original', 'transformed']
    }
    BOOLEAN_PROPERTIES = [
        'nans_to_zeros', 'show_data_range', 'grid', 'use_middle', 'show_number'
    ]
    STRING_PROPERTIES = [
        'file', 'file_type', 'overlay_previous', 'orientation',
        'summary_method', 'title', 'color', 'negative_color', 'transform',
        'y_axis_values', 'type', 'second_file', 'operation', 'link_file',
        'line_style', 'title', 'bw_color', 'link_color'
    ]
    FLOAT_PROPERTIES = {
        'max_value': [-np.inf, np.inf],
        'min_value': [-np.inf, np.inf],
        'log_pseudocount': [-np.inf, np.inf],
        'alpha': [0, 1],
        'height': [0, np.inf],
        'fontsize': [0, np.inf],
        'line_width': [0, np.inf],
        'scale_link_height': [0, np.inf],
        'scale_line_width': [0, np.inf]
    }
    INTEGER_PROPERTIES = {'number_of_bins': [1, np.inf]}

    # The color can only be a color
    # negative_color can only be a color or None

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.bw = pyBigWig.open(self.properties['file'])
        self.bw2 = None
        if 'second_file' in self.properties['operation']:
            if self.properties['second_file'] is None:
                raise InputError(f"operation: {self.properties['operation']}"
                                 " requires to set the parameter"
                                 " second_file.")
            else:
                self.bw2 = pyBigWig.open(self.properties['second_file'])

    def set_properties_defaults(self):
        super(SashimiBigWigTrack, self).set_properties_defaults()
        super(SashimiBigWigTrack, self).process_type_for_coverage_track()
        self.process_color('bw_color')
        self.process_color('link_color')
        self.properties['scale_link_height'] = 1
        self.properties['scale_line_width'] = 1
        if self.properties['negative_color'] is None:
            self.properties['negative_color'] = self.properties['bw_color']
        else:
            self.process_color('negative_color')
        if self.properties['operation'] != 'file':
            self.checkoperation()
            if self.properties['transform'] != 'no':
                raise InputError("'operation' and 'transform' cannot be set at"
                                 " the same time.")
            if self.properties['y_axis_values'] == 'original':
                self.log.warning("*Warning* 'operation' is used and "
                                 "'y_axis_values' was set to 'original'. "
                                 "'y_axis_values' can only be set to "
                                 "'original' when 'transform' is used.\n"
                                 " It will be set as 'transformed'.\n")
                self.properties['y_axis_values'] = 'transformed'

        #FROM LINK
        self.pos_height = None
        self.neg_height = None
        self.interval_tree, min_score, max_score, has_score = self.process_link_file(
            self.properties['region'])
        if self.properties['line_width'] is None and not has_score:
            self.log.warning("*WARNING* for section "
                             f"{self.properties['section_name']}"
                             " no line_width has been set but some "
                             "lines do not have scores."
                             "line_width has been set to "
                             "0.5.\n")
            self.properties['line_width'] = 0.5

        if self.properties['scale_line_width'] is None and has_score:
            self.properties['scale_line_width'] = 2

        self.colormap = None
        # check if the color given is a color map
        is_colormap = self.process_color('link_color',
                                         colormap_possible=True,
                                         default_value_is_colormap=False)
        if is_colormap:
            if not has_score:
                self.log.warning("*WARNING* for section "
                                 f"{self.properties['section_name']}"
                                 " a colormap was chosen but some "
                                 "lines do not have scores."
                                 "Color has been set to "
                                 f"{DEFAULT_LINKS_COLOR}.\n")
                self.properties['link_color'] = DEFAULT_LINKS_COLOR
            else:
                self.colormap = self.properties['link_color']

        if self.colormap is not None:
            if self.properties['min_value'] is not None:
                min_score = self.properties['min_value']
            if self.properties['max_value'] is not None:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score, vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.properties['link_color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def plot(self, ax, chrom_region, start_region, end_region):
        self.pos_height = 0
        self.neg_height = 0
        count = 0

        temp_end_region, temp_nbins, scores_per_bin = self.get_scores(
            'self.bw', self.properties['file'], chrom_region, start_region,
            end_region)
        if scores_per_bin is None:
            self.log.warning(
                "Scores could not be computed. This will generate an empty track\n"
            )
            return

        if self.properties['nans_to_zeros'] and np.any(
                np.isnan(scores_per_bin)):
            scores_per_bin[np.isnan(scores_per_bin)] = 0

        x_values = np.linspace(start_region, temp_end_region, temp_nbins)
        # compute the operation
        operation = self.properties['operation']
        # Substitute log by np.log to make it evaluable:
        operation = operation.replace('log', 'np.log')
        if operation == 'file':
            pass
        elif 'second_file' not in operation:
            try:
                new_scores_per_bin = eval('[' + operation +
                                          ' for file in scores_per_bin]')
                new_scores_per_bin = np.array(new_scores_per_bin)
            except Exception as e:
                raise Exception("The operation in section "
                                f"{self.properties['section_name']} could not "
                                f"be computed: {e}")
            else:
                scores_per_bin = new_scores_per_bin
        else:
            temp_end_region2, temp_nbins2, scores_per_bin2 = self.get_scores(
                'self.bw2', self.properties['second_file'], chrom_region,
                start_region, end_region)
            if scores_per_bin2 is None:
                self.log.warning(
                    "Scores for second_file could not be computed. This will generate an empty track\n"
                )
                return

            if self.properties['nans_to_zeros'] and np.any(
                    np.isnan(scores_per_bin2)):
                scores_per_bin2[np.isnan(scores_per_bin2)] = 0

            x_values2 = np.linspace(start_region, temp_end_region2,
                                    temp_nbins2)
            if not np.all(x_values == x_values2):
                raise Exception(
                    'The two bigwig files are not compatible on this region:'
                    f'{chrom_region}:{start_region}-{end_region}')
            # compute the operation
            try:
                new_scores_per_bin = eval('[' + operation +
                                          ' for file, second_file in'
                                          ' zip(scores_per_bin,'
                                          ' scores_per_bin2)]')
                new_scores_per_bin = np.array(new_scores_per_bin)
            except Exception as e:
                raise Exception("The operation in section "
                                f"{self.properties['section_name']} could not "
                                f"be computed: {e}")
            else:
                scores_per_bin = new_scores_per_bin

        transformed_scores = transform(scores_per_bin,
                                       self.properties['transform'],
                                       self.properties['log_pseudocount'],
                                       self.properties['file'])

        plot_coverage(ax, x_values, transformed_scores, self.plot_type,
                      self.size, self.properties['bw_color'],
                      self.properties['negative_color'],
                      self.properties['alpha'], self.properties['grid'])

        plot_ymin, plot_ymax = ax.get_ylim()
        plot_ymax = eval(f'[{operation} for file in [plot_ymax]]')[0]

        # PLOT LINK
        arcs_in_region = sorted(
            self.interval_tree[chrom_region][start_region:end_region])
        for idx, interval in enumerate(arcs_in_region):
            # skip intervals whose start and end are outside the plotted region
            if interval.begin < start_region and interval.end > end_region:
                continue
            score_start = float(
                self.bw.values(chrom_region, interval.begin,
                               interval.begin + 1)[0])
            score_end = float(
                self.bw.values(chrom_region, interval.end,
                               interval.end + 1)[0])

            if operation == "file":
                pass
            else:
                score_start = eval(
                    f'[{operation} for file in [score_start]]')[0]
                score_end = eval(f'[{operation} for file in [score_end]]')[0]

            if self.properties['line_width'] is not None:
                self.line_width = float(self.properties['line_width'])
            else:
                self.line_width = self.properties['scale_line_width'] * np.log(
                    interval.data[4] + 1) * 1.1

            self.show_number = self.properties['show_number']
            self.plot_bezier(ax, interval, idx, score_start, score_end,
                             plot_ymax)
            count += 1

        # this height might be removed
        self.neg_height *= 1
        self.pos_height *= 1.05
        self.log.debug(f"{count} links plotted")

        if self.properties['min_value'] == None:
            ymin = min(plot_ymin, self.neg_height)
        else:
            ymin = min(plot_ymin, self.properties['min_value'],
                       self.neg_height)

        if self.properties['max_value'] == None:
            ymax = max(plot_ymax, self.pos_height)
        else:
            ymax = max(plot_ymax, self.properties['max_value'],
                       self.pos_height)

        ymax = transform(np.array([ymax]), self.properties['transform'],
                         self.properties['log_pseudocount'], 'ymax')[0]

        ymin = transform(np.array([ymin]), self.properties['transform'],
                         self.properties['log_pseudocount'], 'ymin')[0]

        if self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, ymin)
        else:
            ax.set_ylim(ymin, ymax)

        return ax

    def plot_bezier(self, ax, interval, idx, start_height, end_height, ymax):

        def cubic_bezier(pts, t):
            b_x = (1 - t)**3 * pts[0][0] + 3 * t * (1 - t)**2 * pts[1][
                0] + 3 * t**2 * (1 - t) * pts[2][0] + t**3 * pts[3][0]
            b_y = (1 - t)**3 * pts[0][1] + 3 * t * (1 - t)**2 * pts[1][
                1] + 3 * t**2 * (1 - t) * pts[2][1] + t**3 * pts[3][1]
            return ((b_x, b_y))

        width = (interval.end - interval.begin)

        height = ymax * 0.25 * self.properties['scale_link_height']
        if self.colormap:
            # translate score field
            # into a color
            rgb = self.colormap.to_rgba(interval.data[4])
        else:
            rgb = self.properties['link_color']

        # Plot below x-axis
        if idx % 2 != 0:
            pts = [(interval.begin, 0), (interval.begin, -height),
                   (interval.end, -height), (interval.end, 0)]
            midpt = cubic_bezier(pts, 0.5)
            minpt = min(
                [cubic_bezier(pts, x)[1] for x in np.arange(0, 1, 0.05)])
            if minpt < self.neg_height:
                self.neg_height = minpt

            pp1 = mpatches.PathPatch(Path(
                pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                     fc="none",
                                     ec=self.properties['link_color'],
                                     lw=self.line_width,
                                     ls=self.properties['line_style'])
            ax.add_patch(pp1)
            if self.show_number:
                ax.text(midpt[0],
                        midpt[1],
                        round(interval.data[4], 3),
                        fontsize=self.properties['fontsize'],
                        horizontalalignment='center',
                        verticalalignment='center',
                        bbox=dict(facecolor='white', edgecolor='none', pad=0))
        # Plot above
        else:
            pts = [(interval.begin, start_height),
                   (interval.begin, height + start_height),
                   (interval.end, height + end_height),
                   (interval.end, end_height)]

            midpt = cubic_bezier(pts, 0.5)
            maxpt = max(
                [cubic_bezier(pts, x)[1] for x in np.arange(0, 1, 0.05)])
            if maxpt > self.pos_height:
                self.pos_height = maxpt

            pp1 = mpatches.PathPatch(Path(
                pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                     fc="none",
                                     ec=self.properties['link_color'],
                                     lw=self.line_width,
                                     ls=self.properties['line_style'])
            ax.add_patch(pp1)
            if self.show_number:
                ax.text(midpt[0],
                        midpt[1],
                        round(interval.data[4], 3),
                        fontsize=self.properties['fontsize'],
                        horizontalalignment='center',
                        verticalalignment='center',
                        bbox=dict(facecolor='white', edgecolor='none', pad=0))

    # This y axis does not show the negative part, which is only Sashimi links
    def plot_y_axis(self,
                    ax,
                    plot_axis,
                    transform='no',
                    log_pseudocount=0,
                    y_axis='tranformed',
                    only_at_ticks=False):
        """
        Plot the scale of the y axis with respect to the plot_axis
        Args:
            ax: axis to use to plot the scale
            plot_axis: the reference axis to get the max and min.
            transform: what was the transformation of the data
            log_pseudocount:
            y_axis: 'tranformed' or 'original'
            only_at_ticks: False: only min_max are diplayed
                           True: only ticks values are displayed

        Returns:

        """
        if not self.properties.get('show_data_range', True):
            return

        def value_to_str(value):
            # given a numeric value, returns a
            # string that removes unneeded decimal places
            if value % 1 == 0:
                str_value = str(int(value))
            else:
                str_value = f"{value:.1f}"
            return str_value

        def untransform(value, transform, log_pseudocount):
            # given a numeric value, transform and log_pseudocount
            # return the value before the transformation
            if transform == 'log':
                return np.exp(value) - log_pseudocount
            elif transform == 'log2':
                return np.exp2(value) - log_pseudocount
            elif transform == 'log10':
                return np.power(10, value) - log_pseudocount
            elif transform == 'log1p':
                return np.expm1(value)
            elif transform == '-log':
                return np.exp(-value) - log_pseudocount

        ymin, ymax = plot_axis.get_ylim()
        ymin = 0
        # If the ticks are closer than epsilon from the top or bottom
        # The vertical alignment of label is adjusted
        epsilon = (ymax - ymin) / 100
        # When the ymax and ymin are plotted (when there is no grid)
        # The tick is shifted inside of epsilon_pretty
        # To avoid to have only half of the width of the line plotted
        epsilon_pretty = epsilon

        if only_at_ticks:
            # plot something that looks like this:
            # tick3 ┐
            #       │
            # tick2-|
            #       │
            # tick1 ┘
            if ymin < ymax:
                ticks_values = [
                    t for t in plot_axis.get_yticks()
                    if t <= ymax and t >= ymin
                ]
            else:
                ticks_values = [
                    t for t in plot_axis.get_yticks()
                    if t >= ymax and t <= ymin
                ]
                ticks_values.sort(reverse=True)
            labels_pos = ticks_values
            if transform == 'no' or y_axis == 'transformed':
                ticks_labels = [value_to_str(t) for t in ticks_values]
            else:
                # There is a transformation and we want to display original values
                ticks_labels = [
                    value_to_str(untransform(t, transform, log_pseudocount))
                    for t in ticks_values
                ]
        elif transform == 'no' or y_axis == 'transformed':
            # This is a linear scale
            # plot something that looks like this:
            # ymax ┐
            #      │
            #      │
            # ymin ┘
            # adjust the positions such that the lines are plotted complete
            # and not only half of the width of the line.
            ticks_values = [ymin + epsilon_pretty, ymax - epsilon_pretty]
            labels_pos = [ymin, ymax]
            ticks_labels = [value_to_str(v) for v in [ymin, ymax]]
            if y_axis == 'transformed' and transform != 'no':
                if transform == 'log1p':
                    ymid_str = "log(1 + x)"
                else:
                    if log_pseudocount == 0:
                        ymid_str = f"{transform}(x)"
                    else:
                        ymid_str = f"{transform}({log_pseudocount} + x)"

                ax.text(0, (ymax + ymin) / 2,
                        ymid_str,
                        verticalalignment='center',
                        horizontalalignment='right',
                        wrap=True)
        else:
            # There is a transformation and we want to display original values
            if ymin * ymax < 0:
                ymid = 0
            else:
                ymid = (ymin + ymax) / 2
            # plot something that looks like this:
            # ymax ┐
            #      │
            # ymid-|
            #      │
            # ymin ┘
            ticks_values = [ymin + epsilon_pretty, ymid, ymax - epsilon_pretty]
            labels_pos = [ymin, ymid, ymax]
            ticks_labels = [
                value_to_str(untransform(v, transform, log_pseudocount))
                for v in [ymin, ymid, ymax]
            ]

        # The lower label should be verticalalignment='bottom'
        # if it corresponds to ymin
        i = 0
        if (ymin < ymax and ticks_values[i] <= ymin + epsilon) \
           or (ymin > ymax and ticks_values[i] >= ymin + epsilon):
            v_al = 'bottom'
            adjusted_value = labels_pos[i] - epsilon
        else:
            v_al = 'center'
            adjusted_value = labels_pos[i]
        ax.text(-0.2,
                adjusted_value,
                ticks_labels[i],
                verticalalignment=v_al,
                horizontalalignment='right')
        x_pos = [0, 0.5]
        y_pos = [ticks_values[i]] * 2
        for i in range(1, len(ticks_values) - 1):
            ax.text(-0.2,
                    labels_pos[i],
                    ticks_labels[i],
                    verticalalignment='center',
                    horizontalalignment='right')
            x_pos += [0.5, 0, 0.5]
            y_pos += [ticks_values[i]] * 3

        # The upper label should be verticalalignment='top'
        # if it corresponds to ymax
        i = len(ticks_values) - 1
        if (ymin < ymax and ticks_values[i] >= ymax - epsilon) \
           or (ymin > ymax and ticks_values[i] <= ymax - epsilon):
            v_al = 'top'
        else:
            v_al = 'center'
        ax.text(-0.2,
                labels_pos[i],
                ticks_labels[i],
                verticalalignment=v_al,
                horizontalalignment='right')
        x_pos += [0.5, 0]
        y_pos += [ticks_values[i]] * 2

        # Finally plot the line:
        ax.plot(x_pos, y_pos, color='black', linewidth=1)

        # Set the lims:
        ax.set_ylim(plot_axis.get_ylim())
        ax.set_xlim(0, 1)
        ax.patch.set_visible(False)

    def process_link_file(self, plot_regions):
        # the file format expected is similar to file format of links in
        # circos:
        # chr1 100 200 chr1 250 300 0.5
        # where the last value is a score.

        if plot_regions is None:
            file_to_open = self.properties['link_file']
        else:
            # To be sure we do not miss links we will intersect with bed with
            # only chromosomes used in plot_regions
            plot_regions_adapted = [(chrom, 0, HUGE_NUMBER)
                                    for chrom, __, __ in plot_regions]
            file_to_open = temp_file_from_intersect(
                self.properties['link_file'], plot_regions_adapted)

        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        has_score = True
        max_score = float(0)
        min_score = float(0)
        file_h = opener(file_to_open)
        for line in tqdm(file_h.readlines()):
            line_number += 1
            line = to_string(line)
            if line.startswith('browser') or line.startswith(
                    'track') or line.startswith('#'):
                continue
            try:
                chrom1, start1, end1, chrom2, start2, end2 = line.strip(
                ).split('\t')[:6]
            except Exception as detail:
                raise InputError('File not valid. The format is chrom1'
                                 ' start1, end1, '
                                 f'chrom2, start2, end2\nError: {detail}\n'
                                 f' in line\n {line}')
            if chrom1 != chrom2:
                self.log.warning(
                    f"Only links in same chromosome are used. Skipping line\n{line}\n"
                )
                continue

            try:
                score = line.strip().split('\t')[6]
            except IndexError:
                has_score = False
                score = np.nan

            try:
                start1 = int(start1)
                end1 = int(end1)
                start2 = int(start2)
                end2 = int(end2)
            except ValueError as detail:
                raise InputError(
                    f"Error reading line: {line_number}. One of the fields is not "
                    f"an integer.\nError message: {detail}")

            assert start1 <= end1, f"Error in line #{line_number}, end1 larger than start1 in {line}"
            assert start2 <= end2, f"Error in line #{line_number}, end2 larger than start2 in {line}"

            if has_score:
                try:
                    score = float(score)
                except ValueError as detail:
                    self.log.warning(
                        f"Warning: reading line: {line}. The score is not valid {score} will not be used. "
                        f"\nError message: {detail}\n")
                    score = np.nan
                    has_score = False
                else:
                    if score < min_score:
                        min_score = score
                    if score > max_score:
                        max_score = score

            if chrom1 not in interval_tree:
                interval_tree[chrom1] = IntervalTree()

            if start2 < start1:
                start1, start2 = start2, start1
                end1, end2 = end2, end1

            if self.properties['use_middle']:
                mid1 = (start1 + end1) / 2
                mid2 = (start2 + end2) / 2
                interval_tree[chrom1].add(
                    Interval(mid1, mid2, [start1, end1, start2, end2, score]))
            else:
                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(
                    Interval(start1, end2,
                             [start1, end1, start2, end2, score]))
            valid_intervals += 1

        if valid_intervals == 0:
            self.log.warning(
                f"No valid intervals were found in file {self.properties['link_file']}.\n"
            )

        file_h.close()
        return (interval_tree, min_score, max_score, has_score)

    def get_scores(self, bw_var, file, chrom_region, start_region, end_region):
        bw = eval(bw_var)
        scores_per_bin = None
        if chrom_region not in bw.chroms().keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in bw.chroms().keys():
                self.log.warning("*Warning*\nNeither " + chrom_region_before +
                                 " nor " + chrom_region + " exists as a "
                                 "chromosome name inside the bigwig file. "
                                 "No score will be computed for"
                                 f" {file}.\n")
                scores_per_bin = np.array([np.nan] *
                                          self.properties['number_of_bins'])

        if scores_per_bin is None and start_region > bw.chroms()[chrom_region]:
            self.log.warning(
                "*Warning*\nThe region to plot starts beyond the"
                " chromosome size. No score will be computed for"
                f" {file}.\n"
                f"{chrom_region} size: {bw.chroms()[chrom_region]}"
                f". Region to plot {start_region}-{end_region}\n")
            scores_per_bin = np.array([np.nan] *
                                      self.properties['number_of_bins'])

        if scores_per_bin is None and end_region > bw.chroms()[chrom_region]:
            self.log.warning(
                "*Warning*\nThe region to plot extends beyond the"
                " chromosome size. Please check.\n"
                f"{chrom_region} size: {bw.chroms()[chrom_region]}"
                f". Region to plot {start_region}-{end_region}\n")
            temp_end_region = bw.chroms()[chrom_region]
            temp_nbins = int(self.properties['number_of_bins'] *
                             (temp_end_region - start_region) /
                             (end_region - start_region))
        else:
            temp_end_region = end_region
            temp_nbins = self.properties['number_of_bins']
        # on rare occasions pyBigWig may throw an error, apparently caused by a corruption
        # of the memory. This only occurs when calling trackPlot from different
        # processors. Reloading the file solves the problem.
        if scores_per_bin is None:
            num_tries = 0
            while num_tries < 5:
                num_tries += 1
                try:
                    scores_per_bin = np.array(
                        bw.stats(
                            chrom_region,
                            start_region,
                            temp_end_region,
                            nBins=temp_nbins,
                            type=self.properties['summary_method'])).astype(
                                float)
                except Exception as e:
                    bw = pyBigWig.open(self.properties['file'])

                    self.log.warning("error found while reading bigwig scores "
                                     f"({e}).\nTrying again."
                                     f" Iter num: {num_tries}.\n")
                    pass
                else:
                    if num_tries > 1:
                        self.log.warning(
                            f"After {num_tries} the scores could be computed.\n"
                        )
                    break
        return temp_end_region, temp_nbins, scores_per_bin

    def __del__(self):
        try:
            self.bw.close()
        except AttributeError:
            pass
        if self.bw2 is not None:
            self.bw2.close()
