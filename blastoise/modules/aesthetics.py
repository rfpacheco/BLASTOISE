BOX_WIDTH_DEFAULT = 80


def print_message_box(message: str, width: int = BOX_WIDTH_DEFAULT) -> None:
    """
    Print ``message`` centered inside an ASCII art box.

    Parameters
    ----------
    message : str
        Text to display.
    width : int, optional
        Number of ``=`` characters in the top/bottom borders.
        If smaller than the message length, it is automatically enlarged.
    """
    width = max(width, len(message))  # ensure the message fits
    border = "=" * width

    print()  # leading blank line
    print(f"|{border}|")
    print(f"|{message:^{width}}|")
    print(f"|{border}|")
    print()  # trailing blank line
