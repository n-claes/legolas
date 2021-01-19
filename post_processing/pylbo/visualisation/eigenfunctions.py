from pylbo.utilities.toolbox import add_pickradius_to_item


class EigenfunctionHandler:
    def __init__(self, data):
        self.data = data
        self._selected_idxs = {}

    def on_point_pick(self, event):
        artist = event.artist
        # if artist is a legend item, return (this attribute has been set manually)
        if hasattr(artist, "is_legend_item"):
            return
        fig, ax = artist.figure, artist.axes
        # This retrieves the indices of the clicked points. Multiple indices are
        # possible depending on an overlapping pickradius. Look which point corresponds
        # to the smallest distance to the mouse click.
        idxs = event.ind
        xdata = artist.get_xdata()
        ydata = artist.get_ydata()
        if len(idxs) == 1:
            idx = idxs[0]
        else:
            mouse_x = event.mouseevent.xdata
            mouse_y = event.mouseevent.ydata
            distances = (mouse_x - xdata[idxs])**2 + (mouse_y - ydata[idxs])**2
            idx = idxs[distances.argmin()]
        xdata = xdata[idx]
        ydata = ydata[idx]
        # handle left clicking
        if event.mouseevent.button == 1:
            # skip if point index is alreadt in list
            if str(idx) in self._selected_idxs.keys():
                return
            marked_point, = ax.plot(
                xdata, ydata, "rx", markersize=8, lw=2, label="marked_point",
            )
            add_pickradius_to_item(item=marked_point, pickradius=1)
            self._selected_idxs.update({f"{idx}": marked_point})
        # handle right clicking
        elif event.mouseevent.button == 3:
            # remove selected index from list
            selected_artist = self._selected_idxs.pop(str(idx), None)
            if selected_artist is not None:
                selected_artist.remove()
        fig.canvas.draw()

    def on_key_press(self, event):
        pass
