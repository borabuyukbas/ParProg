#include "plotting/plotter.h"
#include "io/image_parser.h"

#include <exception>

void Plotter::write_and_clear(){
    // create plot serial number string
    std::string serial_number_string = std::to_string(image_serial_number);
    while(serial_number_string.length() < 9){
        serial_number_string = "0" + serial_number_string;
    }

    std::string file_name = filename_prefix + "_" + serial_number_string + ".bmp";
    ImageParser::write_bitmap(output_folder_path / file_name, image);
    clear_image();
    image_serial_number += 1;
}

BitmapImage::BitmapPixel Plotter::get_pixel(std::uint32_t x, std::uint32_t y){
    return image.get_pixel(y, x);
}

void Plotter::highlight_position(Vector2d<double> position, std::uint8_t red, std::uint8_t green, std::uint8_t blue){
    // check for bounding box range
    if (plot_bounding_box.contains(position))
    {
        // // fill in cross across entire image
        uint32_t image_width = image.get_width();
        uint32_t image_height = image.get_height();
       
        double x_min = plot_bounding_box.x_min;
        double x_max = plot_bounding_box.x_max;
        
        double y_min = plot_bounding_box.y_min;
        double y_max = plot_bounding_box.y_max;

        double universe_width = x_max - x_min;
        double universe_height = y_max - y_min;

        // from universe to image space
        double x_pixel = (position[0] - x_min) / universe_width * (image.get_width() - 1);
        double y_pixel = (position[1] - y_min) / universe_height * (image.get_height() - 1);

        // draw horizontal line
        for (int i = 0; i < image_width; i++)
        {
            mark_pixel(i, y_pixel, red, green, blue);
        }

        // draw vertical line
        for (int i = 0; i < image_height; i++)
        {
            mark_pixel(x_pixel, i, red, green, blue);
        }
    }
}