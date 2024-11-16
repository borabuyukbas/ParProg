#pragma once
#include "structures/bounding_box.h"
#include "image/bitmap_image.h"
#include "image/pixel.h"
#include "structures/universe.h"
#include <cstdint>
#include <set>
#include <filesystem>
#include "structures/vector2d.h"

class Plotter{
public:
    Plotter(BoundingBox bb, const std::filesystem::path & arg_output_folder_path, std::uint32_t plot_width_arg, std::uint32_t plot_height_arg) :
        plot_bounding_box(bb), 
        output_folder_path(arg_output_folder_path),
        plot_width(plot_width_arg),
        plot_height(plot_height_arg),
        image(BitmapImage(plot_height_arg, plot_width_arg)),
        image_serial_number(0){
        // set default filename prefix
        filename_prefix = "plot";
    }

    void add_bodies_to_image(Universe& universe);
    void highlight_position(Vector2d<double> position, std::uint8_t red, std::uint8_t green, std::uint8_t blue);
    
    void set_plot_bounding_box(BoundingBox bb){
        plot_bounding_box = bb;
    }

    void write_and_clear();
    
    void clear_image(){
        image = BitmapImage(plot_height, plot_width);
    }

    void set_filename_prefix(std::string prefix){
        filename_prefix = prefix;
    }

    void mark_position(Vector2d<double> position, std::uint8_t red, std::uint8_t green, std::uint8_t blue) {
        if (!plot_bounding_box.contains(position)) return;

        // calculate pixel coordinates using relative width/height
        std::uint32_t x = ((position[0] - plot_bounding_box.x_min) / (plot_bounding_box.x_max - plot_bounding_box.x_min)) * (image.get_width() - 1);
        std::uint32_t y = ((position[1] - plot_bounding_box.y_min) / (plot_bounding_box.y_max - plot_bounding_box.y_min)) * (image.get_height() - 1);
        BitmapImage::BitmapPixel pixel(red, green, blue);

        image.set_pixel(y, x, pixel);
    }
    void mark_pixel(std::uint32_t x, std::uint32_t y, std::uint8_t red, std::uint8_t green, std::uint8_t blue) {
        // exception handling
        if (x >= plot_width || y >= plot_height) {
            throw std::out_of_range("Pixel's position is outside the image bounds!");
    }
        BitmapImage::BitmapPixel pixel(red, green, blue);
        image.set_pixel(y, x, pixel);
    }

    BitmapImage::BitmapPixel get_pixel(std::uint32_t x, std::uint32_t y);

    std::uint32_t get_next_image_serial_number(){
        return image_serial_number;
    }

private:
    std::string filename_prefix;
    std::uint32_t image_serial_number;
    BitmapImage image;
    BoundingBox plot_bounding_box;
    std::uint32_t plot_width, plot_height;
    std::filesystem::path output_folder_path;
};